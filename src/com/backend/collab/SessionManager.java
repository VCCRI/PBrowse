package com.backend.collab;

import java.security.MessageDigest;
import java.security.SecureRandom;
import java.sql.SQLException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Properties;
import java.io.FileInputStream;

import org.eclipse.jetty.websocket.api.Session;

import com.database.mysql.DBA;
import com.database.mysql.DataDesc;
import com.database.mysql.Comment;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonPrimitive;

import static com.database.mysql.DBA.dba;

/**
 * Manages sessions
 * Contains all functions called by WebSocketHandler
 * @author root
 *
 * A singleton object
 */

public class SessionManager 
{
	//global debug flag
	public static final boolean DEBUG = false;
	
	//Set of all active web sockets
	private Set<Session> wss; 
	
	//Set of all active sessions, each session contains >=1 user
	private Set<PSession> sessions;
	
	//map connecting individual users -> sessions
	private HashMap<IUser,PSession> user_session;
	
	//map connecting web socket -> user
	private HashMap<Session,IUser> ws_user;
		
	//map connecting uuid to session -> user
	private HashMap<String,PSession> uuid_session;
		
	//map of logged in users, username -> IUser
	private HashMap<String,IUser> username_user;
	
	//session manager is a  singleton instance
	private static SessionManager ism = null;
	
	//secure number generator used for creating authentication tokens
	private SecureRandom sr = null;
	
	//mailer thread
	private final SendMailTLS send_mail = new SendMailTLS();
	private Thread mail_thread;
	
	//domain specification used in email messages to direct users back to the application
	private String DOMAIN = "";
	
	private SessionManager()
	{
		wss = new HashSet<Session>();
		sessions = new HashSet<PSession>();
		user_session = new HashMap<IUser,PSession>();
		ws_user = new HashMap<Session,IUser>();
		uuid_session = new HashMap<String,PSession>();
		username_user = new HashMap<String,IUser>();
		sr = new SecureRandom();

		//read email configuration
		Properties p = new Properties();
		try {
			p.load(new FileInputStream("/opt/pbrowse/email_bot_config.properties"));
		} catch (Exception e) {
			e.printStackTrace();
		}

		DOMAIN = p.getProperty("domain");
		
		//start the mail thread
		mail_thread = new Thread(send_mail);
		mail_thread.start();
	}
	
	public static SessionManager getISM()
	{
		if (ism == null)
		{
			ism = new SessionManager();
		}
		return ism;
	}
	
	/**
	 * Add a new session, corresponding to a new user
	 * @param s the websocket session
	 */
	public void openSocket(Session s)
	{
		//make a new user, the leader of the new session
		IUser user = new IUser(s,true);

		//map ws session to user
		wss.add(s);
		ws_user.put(s, user);
		user_session.put(user, null);
		
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "new-connection");
		msg.addProperty("status", "connected");
		user.sendMessage(msg.toString());
	}
	
	/**
	 * Sends the message to all active sessions
	 * @param msg - message to send
	 */
	public void broadcastMessage(String msg)
	{
		for (Session s : wss)
		{
			try
			{
				s.getRemote().sendStringByFuture(msg);
			}
			catch (Exception e)
			{
				//the socket is maybe closed
			}
		}
	}
	
	/**
	 * Performs SHA512 hashing on the input string, producing a 128 character hex string
	 * @param s - input string
	 * @return constant length hashed string
	 */
	public String sha512(String s)
	{
		MessageDigest md;
		String hash = null;
		try {
			md = MessageDigest.getInstance("SHA-512");
			md.update(s.getBytes("UTF-8"));
			byte[] digest = md.digest();
			hash = String.format("%064x", new java.math.BigInteger(1, digest));
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return hash;
	}
	
	/**
	 * logs in a user if the provided credentials are valid. Generates a uauth token as well.
	 * Also returns the users active group memberships
	 * @param s - the websocket session
	 * @param params - username and password
	 */
	public void doLogin(Session s, JsonObject params)
	{
		String username = params.get("username").getAsString();
		String password = params.get("password").getAsString();
		IUser u = ws_user.get(s);
		
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "login");
		
		password = this.sha512(password);
		
		HashMap<String, String> props = null;
		try {
			props = dba().doLogin(username, password);
		} catch (SQLException e) {
			e.printStackTrace();
		}
		if (props == null)
		{
			msg.addProperty("status", "fail");
			msg.addProperty("error", "Bad username or password.");
			u.sendMessage(msg.toString());
			return;
		}
		
		if (username_user.containsKey(username))
		{
			//logout previous user - if they still exist.
			IUser ou = username_user.get(username);
			this.doLogout(ou.getSession());
		}
		
		//register username
		username_user.put(username, u);
		
		u.username = username;
		u.nick = props.get("nick");
		u.loggedIn = true;
		
		msg.addProperty("status", "success");
		msg.addProperty("nick", u.nick);
		msg.addProperty("username", u.username);
		
		try {
			msg.addProperty("saved_config", dba().getUserSavedConfig(u.username));
		} catch (SQLException e2) {
			//noop failure
		}
		
		//generate an authentication token for the user's cookie when uploading files
		String pre = "";
		for (int i=0; i<128; i++)
		{
			pre += ""+sr.nextLong();
		}
		pre += username;
		
		String hash = this.sha512(pre);
		
		u.setUAuth(hash);
		msg.addProperty("uauth", hash);
		
		//keep the sessionkey in the database in case the server ever crashes
		try {
			dba().storeSessionKey(username, hash);
		} catch (SQLException e1) {
			//cannot fail for any reason
			e1.printStackTrace();
		}
		//---------------------------------------------------------------------------
		
		//get user group affiliations
		JsonObject g;
		try {
			g = dba().getUserGroups(u.username);
		} catch (SQLException e) {
			e.printStackTrace();
			return;
		}
		msg.add("groups", g);
		u.membership = g;
		
		u.sendMessage(msg.toString());
		return;
	}
	/**
	 * Logs out the calling user
	 * @param s - the websocket session
	 */
	public void doLogout(Session s)
	{
		IUser u = ws_user.get(s);
		if (u == null)
		{
			return;
		}
		
		if (u.loggedIn == false)
		{
			return;
		}
		username_user.remove(u.username);
		
		//instructs client to discard cookies
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "logout");
		u.sendMessage(msg.toString());
		
		//force user to leave session
		if (user_session.get(u) != null)
			this.leaveSession(s);

		//cleans all login related data from user
		u.clearData();
	}
	/**
	 * Registers a new account with the given parameters. The account is initially inactive and must
	 * be activated by the user though the email verification process before it can be logged in to.
	 * A auth token is generated and sent as an email to the user's email address. Registration fails
	 * if the username is already taken.
	 * @param s - the websocket session
	 * @param params - username, password, nick, email
	 */
	public void doRegister(Session s, JsonObject params)
	{
		IUser u = ws_user.get(s);
		String username, password, nick, email;
		
		username = params.get("username").getAsString();
		password = params.get("password").getAsString();
		email = params.get("email").getAsString();
		nick = params.get("nick").getAsString();
		
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "register");
		
		password = this.sha512(password);
		params.addProperty("password", password);
		
		try {
			if (dba().doRegister(username, nick, email, password))
			{
				msg.addProperty("status", "success");
				u.sendMessage(msg.toString());
				
				String token = this.generateToken(u.username).substring(0,32);
				dba().setToken(username, token, DBA.FUNCTION_REGISTERACCOUNT, params.toString());
				
				send_mail.enqueueMessage("New Account", email, 
						"Thank you for registering an account with PBrowse, the Collaborative genome browser<br>"
						+ "Follow this link to activate your account:<br><br>"
						+ DOMAIN+"?token="+token+"<br><br>"
						+ "Regards,<br>"
						+ "PBrowse Team.");
			}
			else
			{
				msg.addProperty("status", "fail");
				msg.addProperty("error", "Username is already taken.");
				u.sendMessage(msg.toString());
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Gets a list of all the user's uploaded files and sends it to the calling user. 
	 * @param s - the websocket session
	 */
	public void getUserFiles(Session s)
	{
		IUser u = ws_user.get(s);
		
		JsonArray ddlist = new JsonArray();
		
		PSession ps = user_session.get(u);
		//if the user is not in a session
		if (ps == null)
		{
			try {
				for (DataDesc dd : dba().getDataFiles(u.username))
				{
					ddlist.add(dd.getJson());
				}
			} catch (SQLException e) {
				e.printStackTrace();
			}
		}
		else
		{
			for (String user : user_session.get(u).getOwners())
			{
				try {
					for (DataDesc dd : dba().getDataFiles(user))
					{
						ddlist.add(dd.getJson());
					}
				} catch (SQLException e) {
					e.printStackTrace();
				}
			}
		}
		
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "user-files");
		msg.add("files", ddlist);
		u.sendMessage(msg.toString());
	}
	
	/**
	 * Returns the IUser instance for a valid logged in user, via their username
	 * @param username
	 * @return
	 */
	public IUser getUserFromName(String username)
	{
		if (username_user.containsKey(username))
		{
			return username_user.get(username);
		}
		return null;
	}
	/**
	 * Gets the interactive session of a given user
	 * @param u - the IUser instance
	 * @return the parallel session instance
	 */
	public PSession getUserSession(IUser u)
	{
		if (user_session.containsKey(u))
		{
			return user_session.get(u);
		}
		return null;
	}
	
	/**
	 * Generate a random session id, not guaranteed unique
	 * @return session id as string
	 */
	private String generateSID(int n)
	{
		String result = "";
		for (int i=0; i<n; i++)
		{
			Character c = (char) (sr.nextInt(26)+65);
			result += c;
		}
		
		return result;
	}
	
	/**
	 * Gets a list of publicly visible sessions and sends it to the calling user
	 * @param s - the websocket session
	 */
	public void getPublicSessions(Session s)
	{
		IUser u = ws_user.get(s);

		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "public-sessions");

		JsonArray list = new JsonArray();
		for (PSession ps : sessions)
		{
			if (ps.getOption("private").equals("true"))
			{	
				continue;
			}
			JsonObject entry = new JsonObject();
			entry.addProperty("uuid", ps.getSID());
			entry.addProperty("numusers", ps.numUsers());
			entry.addProperty("name", ps.getOption("name"));
			
			boolean hascode = ps.getOption("code").length() != 0;
			entry.addProperty("code", hascode);
			
			list.add(entry);
		}
		msg.add("sessions", list);
		
		u.sendMessage(msg.toString());
	}

	/**
	 * make a new interactive session with the calling user as leader. The session is
	 * identified by its randomly generated unique uuid. Call fails if the user is already
	 * part of a session
	 * @param s - the websocket session
	 * @param params - private
	 * @return id of new session, or old id if session exists
	 */
	public String newSession(Session s, JsonObject params)
	{
		Boolean isprivate = false;
		isprivate = params.get("private").getAsBoolean();
		if (isprivate == null)
		{
			System.out.println("Param missing.");
			return null;
		}
		
		IUser u = ws_user.get(s);
		String uuid = null;
		
		if (user_session.get(u) != null)
		{
			uuid = user_session.get(u).getUUID();
		}
		else
		{
			PSession is = new PSession(u);

			String password = params.get("code").getAsString();
			if (password != null)
			{
				is.setOption("code", password);
			}
			is.setOption("private", isprivate?"true":"false");
			String sessionname = params.get("name").getAsString();
			if (sessionname != null)
			{
				is.setOption("name", sessionname);
			}
			
			uuid = generateSID(6);
			while (uuid_session.containsKey(uuid) == true)
			{
				uuid = generateSID(6);
			}
            is.setSID(uuid);
            
            is.genome = params.get("genome").getAsString();
            
            if (is.genome.equals("generic"))
            {
            	is.custom_genome = params.get("custom_genome").getAsString();
            }
			
			//store the new session
			sessions.add(is);
			//map user to interactive session
			user_session.put(u, is);
			//map uuid to session
			uuid_session.put(uuid, is);
			
			//set tiers of session
			is.setTiers(params.getAsJsonArray("tiers"));
			u.updateUserPos(params.getAsJsonObject("status"));
			
			System.out.println("New session uuid:"+uuid+" code:"+password+" isprivate:"+isprivate);
		}
		
		//respond to user with uuid of interactive session
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "new-session");
		msg.addProperty("session", uuid);
		u.sendMessage(msg.toString());
		
		return uuid;
	}
	
	/**
	 * Removes session, all users are kicked from session
	 * @param s - the websocket session
	 * @return success if session closed, false otherwise
	 */
	public boolean closeSession(Session s)
	{
		IUser u = ws_user.get(s);
		if (user_session.get(u) == null || u.isLeader == false)
		{
			return false;
		}
		
		PSession is = user_session.get(u);
		
		//acknowledge
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "close-session");
		is.updateAllUsers(msg.toString());

		//end session for 'all' users
		for (IUser su: is.userlist())
		{
			user_session.put(su, null);
		}
		//force close
		is.closeSession();
		uuid_session.put(is.getUUID(), null);
		
		sessions.remove(is);
		u.isLeader = false;
		return true;
	}
	
	/**
	 * The calling user joins an interactive session identified by the uuid parameter. Sends
	 * a message to the user updating them on the immediate status of all users already
	 * in the session. Call fails if the user is already part of a session
	 * @param s - the websocket session
	 * @param params - uuid
	 * @return true if the session was joined
	 */
	public boolean joinSession(Session s, JsonObject params)
	{
		IUser u = ws_user.get(s);
		if (!u.loggedIn)
		{
			return false;
		}
		
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "join-session");
		PSession is = null;
		boolean success = true;
		
		String uuid = null;
		uuid = params.get("uuid").getAsString();
		if (uuid == null)
		{
			success = false;
			msg.addProperty("error","No session id provided.");
		}
		else
		{
			is = uuid_session.get(uuid);
			if (is == null)
			{
				success = false;
				msg.addProperty("error","Missing entry code or non-existent session.");
			}
			else
			{
				String code = params.get("code").getAsString();
				if (code != null && is.getOption("code") != null)
				{
					if (!code.equals(is.getOption("code")))
					{
						success = false;
						msg.addProperty("error","Missing entry code or non-existent session.");
					}
				}
			}
		}
		if (user_session.containsKey(u))
		{
			if (user_session.get(u) != null)
			{
				//user has to leave a session before joining a new one
				msg.addProperty("error","User belongs to session already.");
				success = false;
			}
		}

		if (!success)
		{
			msg.addProperty("status", "fail");
			u.sendMessage(msg.toString());
			return false;
		}
		
		//update their profile with session-join details
		u.updateUserPos(params);
		
		//enter user into the session
		if (!is.addUser(u))
		{
			msg.addProperty("status", "fail");
			msg.addProperty("error","User has been blacklisted from this session!");
			u.sendMessage(msg.toString());
			return false;
		}
		user_session.put(u, is);
		
		//respond to user with uuid of interactive session
		msg.addProperty("status", "connected");
		msg.addProperty("session", uuid);
		msg.addProperty("joiner", u.username);
		msg.addProperty("genome", is.genome);
		msg.addProperty("custom_genome", is.custom_genome);
		msg.add("options", is.getOptions());
		msg.add("tiers", is.getTiers());
		msg.add("user-status", is.getUsersStatus());
		msg.add("chathistory",is.chatHistory);
		msg.add("blacklist", is.getBlacklist());
		
		//get files of all users
		JsonArray ddlist = new JsonArray();
		for (String user : user_session.get(u).getOwners())
		{
			try {
				for (DataDesc dd : dba().getDataFiles(user))
				{
					ddlist.add(dd.getJson());
				}
			} catch (SQLException e) {
				System.out.println("Could not retrieve files of: "+user);
			}
		}
		msg.add("allfiles", ddlist);
		
		is.updateAllUsers(msg.toString());
		return true;
	}
	
	/**
	 * calling user leaves the interactive session. Call fails if the user
	 * is not part of a session
	 * @param s - the websocket session
	 * @return success if the user was removed from the session
	 */
	public boolean leaveSession(Session s)
	{
		IUser u = ws_user.get(s);
		
		if (!user_session.containsKey(u))
		{
			//user not part of a session
			System.out.println("User not part of a session");
			return false;
		}
		
		PSession is = user_session.get(u);
		user_session.remove(u);
		
		//acknowledge
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "leave-session");
		msg.addProperty("username", u.username);
		msg.addProperty("status", "disconnected");
		
		//the user may have disconnected
		u.sendMessage(msg.toString());
		
		if (is != null)
		{
			//if removing the user empties the session, delete the session
			if (is.removeUser(u) == 0)
			{
				sessions.remove(is);
				uuid_session.put(is.getUUID(), null);
			}
			else
			{
				msg.add("user-status", is.getUsersStatus());
				is.updateAllUsers(msg.toString());
			}
		}
		return true;
	}
	
	/**
	 * The websocket connection is being closed, the user is removed from their shared session
	 * if any
	 * @param s - the websocket session
	 */
	public void closeSocket(Session s)
	{
		//leave session per normal
		IUser u = ws_user.get(s);

		if (u != null)
		{
			if (user_session.get(u) != null)
				this.leaveSession(s);
		}
		
		ws_user.remove(s);
		wss.remove(s);
		u.setSession(null);
	}
	
	/**
	 * Attempt to resume the previously logged in session of a calling user. The client stores
	 * this information as part of their cookies, and sends it to the server the first time they
	 * open a connection. If the credentials are valid, the user is immediately logged back
	 * into the system.
	 * @param s - the websocket session
	 * @param params - username, uauth token
	 */
	public void doResumeSession(Session s, JsonObject params)
	{
		params.addProperty("rtype", "login");

		String username = params.get("username").getAsString();
		String uauth = params.get("uauth").getAsString();
		
		//the user we are resuming
		IUser u = username_user.get(username);
		if (u != null)
		{
			if (u.getUAuth().equals(uauth))
			{
				//valid data - connect this new session to the user
				u.setSession(s);
				ws_user.put(s, u);
				
				//respond with login type message
				params.addProperty("status", "success");
				params.addProperty("nick", u.nick);
				
				JsonObject g = null;
				try {
					g = dba().getUserGroups(u.username);
					params.addProperty("saved_config", dba().getUserSavedConfig(u.username));
				} catch (SQLException e) {
					e.printStackTrace();
					return;
				}
				params.add("groups", g);
			}
		}
		else
		{
			//the server cannot recall the user - might be due to reset or anything else
			//attempt to check the sessionkey from the database
			
			JsonObject result = null;
			try {
				result = dba().checkSessionKey(username, uauth);
			} catch (SQLException e) {
				e.printStackTrace();
			}

			if (result != null)
			{
				IUser caller = ws_user.get(s);
				
				//register user
				username_user.put(username, caller);
				caller.loggedIn = true;
				caller.username = result.get("username").getAsString();
				caller.nick = result.get("nick").getAsString();
				caller.setUAuth(uauth);
				
				//get user group affiliations
				JsonObject g = null;
				try {
					g = dba().getUserGroups(caller.username);
					params.addProperty("saved_config", dba().getUserSavedConfig(caller.username));
				} catch (SQLException e) {
					e.printStackTrace();
					return;
				}
				caller.membership = g;
				
				//create response
				params.add("groups", g);
				params.addProperty("status", "success");
				params.addProperty("username", caller.username);
				params.addProperty("nick", caller.nick);
			}
			else
			{
				params.addProperty("status", "resumefail");
				params.addProperty("error", "Failed to resume session. Bad session data.");
			}
		}
		
		//since we can't guarantee the user will exist, send message via the provided session
		try {
			s.getRemote().sendStringByFuture(params.toString());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Sends a message to all the users of an interactive session
	 * @param s - the websocket session
	 * @param jsonStatus - the set of parameters to forward to all group members
	 */
	public void notifySessionUsers(Session s, JsonObject jsonStatus)
	{
		IUser user = ws_user.get(s);
		PSession is = user_session.get(user);
		
		if (is == null)
		{
			return;
		}

		SimpleDateFormat sdfDate = new SimpleDateFormat("HH:mm:ss");
	    Date now = new Date();
	    String strDate = sdfDate.format(now);
	    jsonStatus.addProperty("time", strDate);
		
		//if this is a position update
		if (jsonStatus.get("chr") != null)
		{
			user.updateUserPos(jsonStatus);
		}
		//session chat message
		else if (jsonStatus.get("smsg") != null)
		{
			jsonStatus.addProperty("username", user.username);
			boolean raw = false;
			if (jsonStatus.get("raw") != null)
			{
				raw = jsonStatus.get("raw").getAsBoolean();
			}
			
			String message;
			if (raw)
			{
				message = jsonStatus.get("smsg").getAsString();
			}
			else
			{
				message = user.username+"("+user.nick+"): "+jsonStatus.get("smsg").getAsString();
			}
			message = "["+strDate+"] - "+message;
			
			JsonPrimitive msg = new JsonPrimitive(message);
			is.chatHistory.add(msg);
		}
		//this is a tier update
		else if (jsonStatus.get("tiers") != null)
		{
			is.setTiers(jsonStatus.getAsJsonArray("tiers"));
			is.genome = jsonStatus.get("genome").getAsString();
			if (is.genome.equals("generic"))
			{
				is.custom_genome = jsonStatus.get("custom_genome").getAsString();
			}
		}
		
		jsonStatus.addProperty("nickname", user.nick);
		jsonStatus.addProperty("caller", user.username);
		//add update rtype to broadcast
		jsonStatus.addProperty("rtype", "update-session-users");
		//determines whether the update invoking user is a leader or not
		jsonStatus.addProperty("isLeader", user.isLeader);
		
		is.updateAllUsers(jsonStatus.toString());
	}

	/**
	 * Update the interactive session options and broadcast the changes to all participants
	 * @param s - the websocket session
	 * @param opt - the options to set, delimited by ':'
	 */
	public void sessionSetOption(Session s, JsonObject opt)
	{
		IUser user = ws_user.get(s);
		if (user.isLeader == false)
		{
			//change nothing, not leader
			return;
		}
		PSession is = user_session.get(user);
		
		//update all options
		String opts[] = opt.get("setoption").getAsString().split(";");
		for (String o : opts)
		{
			String name = o.split(":")[0];
			String val = o.split(":")[1];

			is.setOption(name,val);
		}

		//broadcast changes to all session users
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "update-session-options");
		msg.add("options", is.getOptions());
		is.updateAllUsers(msg.toString());
	}

	/**
	 * Deletes the file indicated by the ID parameter. Call fails if the user does not
	 * own the file - i.e. they were not the one who uploaded it
	 * @param s - the websocket session
	 * @param params - id of the file
	 */
	public void deleteUserFile(Session s, JsonObject params) 
	{
		IUser user = ws_user.get(s);
		
		int id = params.get("id").getAsInt();
				
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "file-delete");
		try {
			if (dba().deleteDataFile(user.username, id))
			{
				msg.addProperty("status", "success");
				msg.addProperty("id", id);
				this.broadcastMessage(msg.toString());
				return;
			}
			else
			{
				msg.addProperty("status", "fail");
				msg.addProperty("error", "File doesn't exist or user does not own it.");
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
		user.sendMessage(msg.toString());
	}

	/**
	 * Toggles the state of a given file between either public or private. Call fails
	 * if the user does not own the file.
	 * @param s - the websocket session
	 * @param params - id of the file
	 */
	public void setUserFilePublicStatus(Session s, JsonObject params) 
	{
		IUser user = ws_user.get(s);
		
		int id = params.get("id").getAsInt();
		int status = params.get("privacy").getAsInt();
				
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "file-public-change");
		try {
			if (dba().setUserFilePublicStatus(user.username, id, status))
			{
				msg.addProperty("status", "success");
				msg.add("dd", dba().getDataFile(id).getJson());
				
				PSession ps = this.getUserSession(user);
				if (ps != null)
					msg.addProperty("uuid", ps.getUUID());
				
				this.broadcastMessage(msg.toString());
			}
			else
			{
				msg.addProperty("status", "fail");
				msg.addProperty("error", "File doesn't exist or user does not own it.");
				user.sendMessage(msg.toString());
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Gets a list of all the publicly accessible files and sends it to the user.
	 * @param s - the websocket session
	 */
	public void getPublicFiles(Session s)
	{
		IUser u = ws_user.get(s);
		
		JsonArray ddlist = new JsonArray();
		try {
			for (DataDesc dd : dba().getPublicFiles())
			{
				ddlist.add(dd.getJson());
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "public-files");
		msg.add("files", ddlist);
		u.sendMessage(msg.toString());
	}
	
	/**
	 * Registers a new comment on a given file, with set parameters. The creation event
	 * is broadcast to all online users. If they are viewing the affected file, they will be
	 * notified of a new comment.
	 * @param s - the websocket session
	 * @param params - trackref (the fileid), chr the chromosome, start position,
	 * end position, ctext (comment text), ispublic (public readability status of the comment)
	 */
	public void makeComment(Session s, JsonObject params)
	{
		IUser user = ws_user.get(s);
		
		long trackref = params.get("trackref").getAsInt();
		String chr = params.get("chr").getAsString();
		long start = params.get("start").getAsInt();
		long end = params.get("end").getAsInt();
		String ctext = params.get("ctext").getAsString();
		int ispublic = params.get("ispublic").getAsInt();

		//create comment object, id is ignored and auto-assigned
		Comment c = new Comment(0,trackref,ctext,chr,start,end,user.username,ispublic);
		
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "make-comment");
		
		int id = 0;
		try {
			id = dba().newComment(c);
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
		if (id != 0)
		{
			//set the id and return it
			c.id = id;
			
			msg.addProperty("status", "success");
			msg.add("comment", c.getJson());
		}
		else
		{
			msg.addProperty("status", "fail");
			msg.addProperty("error", "Stacktrace tells all.");
		}

		
		if (ispublic == DBA.PRIVACY_ONLYME)
		{
			user.sendMessage(msg.toString());
		}
		else if (ispublic == DBA.PRIVACY_PRIVATE)
		{
			PSession ps = this.getUserSession(user);
			if (ps != null)
				ps.updateAllUsers(msg.toString());
			else
				user.sendMessage(msg.toString());
		}
		else
		{
			//update all users about the comment change, they can choose
			//to ignore or act upon the contents if they are viewing the same file
			broadcastMessage(msg.toString());
		}
	}
	
	/**
	 * Retrieves a list of all the comments made on a specified file and sends it to the
	 * calling user.
	 * @param s - the websocket session
	 * @param params - trackref (the file id)
	 */
	public void getAllComments(Session s, JsonObject params)
	{
		IUser user = ws_user.get(s);
		
		JsonElement arr = params.get("trackref");
		if (arr == null)
		{
			return;
		}
		JsonArray tracks = arr.getAsJsonArray();

		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "all-track-comments");
		
		JsonArray cms = new JsonArray();
		
		//get comments for all requested tracks
		for (int i=0; i<tracks.size(); i++)
		{
			long trackref = tracks.get(i).getAsLong();
			try {
				
				//user is part of collaborative session, get all users' comments
				ArrayList<Comment> comments = null;
				PSession ps = user_session.get(user);
				if (ps != null)
				{
					comments = dba().getAllComments(ps.getOwners(),trackref,user.username);
				}
				else 
				{
					List<String> us = new ArrayList<String>(1);
					//call just for this user - if logged in
					if (user.loggedIn)
					{
						us.add(user.username);
						comments = dba().getAllComments(us,trackref,user.username);
					}
					else
					{
						us.add("");
						comments = dba().getAllComments(us,trackref,user.username);
					}
				}
				
				//get json array of data
				for (Comment c : comments)
				{
					cms.add(c.getJson());
				}
			} catch (SQLException e) {
				e.printStackTrace();
			}
		}

		msg.addProperty("status", "success");
		msg.add("comments", cms);
		
		user.sendMessage(msg.toString());
	}
	
	/**
	 * Deletes a comment by its id, the deletion, if successful, is broadcast to all active users.
	 * If they can see the affected comment, it is removed from their view.
	 * @param s - the websocket session
	 * @param params - id (of the comment)
	 */
	public void deleteComment(Session s, JsonObject params)
	{
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		long id = params.get("id").getAsLong();
		long trackref = params.get("trackref").getAsLong();
		
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "delete-comment");
		
		try 	
		{
			dba().deleteComment(id, trackref, user.username);
			msg.addProperty("status", "success");
			msg.addProperty("id", id);
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
			msg.addProperty("status", "fail");
		}

		//update all users about the comment change, they can choose
		//to ignore or act upon the contents
		broadcastMessage(msg.toString());
	}

	/**
	 * Creates a new group with the given name. The calling user is automatically made
	 * owner of the new group.
	 * @param s - the websocket session
	 * @param params - groupname
	 */
	public void createGroup(Session s, JsonObject params) 
	{
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		if (params.get("groupname") == null)
			return;
		
		String groupname = params.get("groupname").getAsString();
		
		boolean success = false;
		try {
			success = dba().createUserGroup(user.username, groupname);
		} catch (SQLException e) {
			//ignore - only 1 possible outcome
		}
		
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "create-group");
		
		if (success)
		{
			msg.addProperty("status", "success");
			msg.addProperty("groupname", groupname);
			
			JsonObject newgroup = new JsonObject();
			newgroup.addProperty("access_level", 1000);
			newgroup.addProperty("count", 1);
			user.membership.add(groupname, newgroup);
		}
		else
		{
			msg.addProperty("status", "fail");
		}
		user.sendMessage(msg.toString());
	}

	/**
	 * Retrieves a list of all the groups the calling user is a member of, along with their
	 * access level for each. Sends this to the user.
	 * @param s - the websocket session
	 */
	public void getUserGroups(Session s) 
	{
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		JsonObject r = null;
		try {
			r = dba().getUserGroups(user.username);
		} catch (SQLException e) {
			e.printStackTrace();
			return;
		}
		
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "user-groups");
		msg.add("groups", r);
		user.sendMessage(msg.toString());
	}

	/**
	 * Gets a list of all the members and files of a particular group and sends it to the calling user
	 * @param s - the websocket session
	 * @param params - groupname
	 */
	public void getGroupInfo(Session s, JsonObject params) 
	{
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		if (params.get("groupname") == null)
			return;
		
		String groupname = params.get("groupname").getAsString();
		
		JsonArray users = null, files = null;
		try {
			users = dba().getGroupUsers(groupname);
			files = dba().getGroupFileDDs(groupname);
		} catch (SQLException e) {
			e.printStackTrace();
			return;
		}
		
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "group-info");
		msg.addProperty("groupname", groupname);
		msg.add("users", users);
		msg.add("files", files);
		user.sendMessage(msg.toString());
	}

	/**
	 * Shares access to a given file with the specified group. The sharer must be the owner
	 * of the file. 
	 * @param s - the websocket session
	 * @param params - groupname, fileid
	 */
	public void shareGroupFile(Session s, JsonObject params) 
	{
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		if (params.get("groupname") == null || params.get("fileid") == null)
		{
			return;
		}
		
		String groupname = params.get("groupname").getAsString();
		long fileid = params.get("fileid").getAsLong();
		
		boolean result = false;
		if (user.membership.get(groupname) != null)
		{
			int user_acl = user.membership
					.getAsJsonObject(groupname)
					.get("access_level").getAsInt();
			
			JsonArray users = null;
			try {
				result = dba().addFileACLToGroup(user.username, fileid, groupname, user_acl);
				users = dba().getGroupUsers(groupname);
			} catch (SQLException e) {
				e.printStackTrace();
				return;
			}
			
			JsonObject msg = new JsonObject();
			msg.addProperty("rtype", "share-group-file");
			
			if (result)
			{
				msg.addProperty("status", "success");
				msg.addProperty("caller", user.username);
				
				//retrieve file DD
				DataDesc dd = null;
				try {
					dd = dba().getDataFile(fileid);
				} catch (SQLException e) {
					e.printStackTrace();
				}
				dd.access_required = user_acl;
				
				msg.add("file", dd.getJson());
				msg.addProperty("groupname", groupname);
			}
			else msg.addProperty("status", "fail");
			
			//notify all group members of new file
			for (JsonElement e : users)
			{
				String name = e.getAsJsonObject()
						.get("username")
						.getAsString();

				System.out.println(name);
				
				if (!username_user.containsKey(name))
				{
					continue;
				}
				
				IUser mem = username_user.get(name);
				mem.sendMessage(msg.toString());
			}
		}
	}

	/**
	 * Invites the specified user into the group and grants them the set of permissions. Notifies
	 * both the caller and the affected user - if they are online, of operation success.
	 * @param s - the websocket session
	 * @param params - groupname, username (the invited user), access_level (new user permissions)
	 */
	public void addUserToGroup(Session s, JsonObject params) 
	{
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		if (params.get("groupname") == null || params.get("username") == null || params.get("access_level") == null)
			return;
		
		String groupname = params.get("groupname").getAsString();
		String username = params.get("username").getAsString();
		int access_level = params.get("access_level").getAsInt();
		
		boolean result = false;
		try {
			result = dba().addUserToGroup(user.username, username, groupname, access_level);
		} catch (SQLException e) {
			e.printStackTrace();
			return;
		}
		
		
		if (result)
		{
			params.addProperty("status", "success");
			
			//tell the added user if they are active
			if (username_user.containsKey(username))
			{
				try {
					int counts = dba().getGroupUserCounts().get(groupname);
					
					//this can be null if the group has no files
					int f_counts = 0;
					HashMap<String,Integer> r = dba().getGroupFileCounts();
					if (r.containsKey(groupname))
						f_counts = r.get(groupname);
					
					params.addProperty("usercount", counts);
					params.addProperty("filecount", f_counts);
				} catch (SQLException e) {
					System.out.println("Failed to retrieve user and file counts.");
				}
				
				params.addProperty("rtype", "added-to-group");
				username_user.get(username).sendMessage(params.toString());
			}
		}
		else params.addProperty("status", "fail");
		
		params.addProperty("rtype", "add-group-user");
		user.sendMessage(params.toString());
	}

	/**
	 * Completely deletes the specified group and removes all associated memberships and file
	 * access lists. Only the group owner can successfully invoke this function. All members of 
	 * the group are notified of its destruction if they are online.
	 * @param s - the websocket session
	 * @param params - groupname
	 */
	public void deleteGroup(Session s, JsonObject params) 
	{
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		if (params.get("groupname") == null)
			return;
		
		String groupname = params.get("groupname").getAsString();
		
		//GET group members FIRST... gawd
		JsonArray users = null;
		boolean result = false;
		try {
			users = dba().getGroupUsers(groupname);
			result = dba().deleteUserGroup(user.username, groupname);
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
		params.addProperty("rtype", "group-deleted");
		
		if (result) 
		{
			params.addProperty("status", "success");
			params.addProperty("caller", user.username);
			
			for (JsonElement e : users)
			{
				String name = e.getAsJsonObject()
						.get("username")
						.getAsString();

				System.out.println(name);
				
				if (!username_user.containsKey(name))
				{
					continue;
				}
				
				IUser mem = username_user.get(name);
				mem.sendMessage(params.toString());
				
				//clear the active membership from all members
				mem.membership.remove(groupname);
			}
		}
		else 
		{
			//otherwise only tell calling user
			params.addProperty("status", "fail");
			user.sendMessage(params.toString());
		}
	}

	/**
	 * Changes the permission of a particular user in the context of a particular group. Notifies 
	 * all group members of any changes.
	 * @param s - the websocket session
	 * @param params - groupname, username, acl (the new permission set)
	 */
	public void updateGroupUserACL(Session s, JsonObject params) 
	{
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		if (params.get("groupname") == null || params.get("username") == null || params.get("acl") == null)
			return;
		
		String groupname = params.get("groupname").getAsString();
		String username = params.get("username").getAsString();
		int acl = params.get("acl").getAsInt();
		
		boolean success = false;
		JsonArray users = null;
		try {
			users = dba().getGroupUsers(groupname);
			success = dba().updateUserAccessLevel(user.username, username, groupname, acl);
		} catch (SQLException e) {
			e.printStackTrace();
			return;
		}
		
		params.addProperty("rtype", "update-group-user-acl");
		
		if (success) 
		{
			params.addProperty("status", "success");
			params.addProperty("caller", user.username);
			
			for (JsonElement e : users)
			{
				String name = e.getAsJsonObject()
						.get("username")
						.getAsString();

				System.out.println(name);
				
				if (!username_user.containsKey(name))
				{
					continue;
				}
				
				IUser mem = username_user.get(name);
				mem.sendMessage(params.toString());
				
				//update membership of affected user only
				if (mem.username.equals(username))
				{
					JsonObject j = mem.membership.getAsJsonObject(groupname);
					j.addProperty("access_level", acl);
					mem.membership.add(groupname, j);
				}
			}
		}
		else 
		{
			//otherwise only tell calling user
			params.addProperty("status", "fail");
			user.sendMessage(params.toString());
		}
	}

	/**
	 * Removes a user from a particular group, notifying all group members. If the username
	 * parameter is not specified, the calling user is removed from the group - equivalent to
	 * leaving. A user is always able to leave a group, unless they are the owner. Calling user
	 * must have user management or ownership privilege to make this change. 
	 * @param s - the websocket session. 
	 * @param params - groupname, username (defaults to self)
	 */
	public void removeGroupUser(Session s, JsonObject params) 
	{
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		if (params.get("groupname") == null)
			return;
		
		String groupname = params.get("groupname").getAsString();

		if (params.get("username") == null)
		{
			//default to calling user
			params.addProperty("username", user.username);
		}
		String username = params.get("username").getAsString();
		
		boolean success = false;
		JsonArray users = null;
		try {
			users = dba().getGroupUsers(groupname);
			success = dba().removeUserFromGroup(user.username, username, groupname);
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
		params.addProperty("rtype", "remove-group-user");
		
		if (success) 
		{
			params.addProperty("status", "success");
			params.addProperty("caller", user.username);
			
			for (JsonElement e : users)
			{
				String name = e.getAsJsonObject()
						.get("username")
						.getAsString();

				System.out.println(name);
				
				if (!username_user.containsKey(name))
				{
					continue;
				}
				
				IUser mem = username_user.get(name);
				mem.sendMessage(params.toString());
				
				//update membership of affected user only
				if (mem.username.equals(username))
				{
					mem.membership.remove(groupname);
				}
			}
		}
		else 
		{
			//otherwise only tell calling user
			params.addProperty("status", "fail");
			user.sendMessage(params.toString());
		}
	}

	/**
	 * Completely removes access of a file from a particular group, notifies all 
	 * the group members. Calling user must have file management or ownership privilege
	 * to make this change.
	 * @param s - the websocket session
	 * @param params - groupname, fileid
	 */
	public void removeGroupFileACL(Session s, JsonObject params) 
	{
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		if (params.get("groupname") == null || params.get("fileid") == null)
			return;
		
		String groupname = params.get("groupname").getAsString();
		long fileid = params.get("fileid").getAsLong();
		
		boolean success = false;
		JsonArray users = null;
		try {
			users = dba().getGroupUsers(groupname);
			success = dba().deleteFileACLFromGroup(user.username, fileid, groupname);
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
		params.addProperty("rtype", "remove-group-file-acl");
		
		if (success) 
		{
			params.addProperty("status", "success");
			params.addProperty("caller", user.username);
			
			for (JsonElement e : users)
			{
				String name = e.getAsJsonObject()
						.get("username")
						.getAsString();

				System.out.println(name);
				
				if (!username_user.containsKey(name))
				{
					continue;
				}
				
				IUser mem = username_user.get(name);
				mem.sendMessage(params.toString());
			}
		}
		else 
		{
			//otherwise only tell calling user
			params.addProperty("status", "fail");
			user.sendMessage(params.toString());
		}
	}
	
	/**
	 * Sends an email message with predetermined subject and user specified content to all
	 * members of a group. Only the group owner may invoke this.
	 * @param s - the websocket session
	 * @param params - groupname, message (the message text)
	 */
	public void sendGroupEmail(Session s, JsonObject params)
	{
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		if (params.get("groupname") == null || params.get("message") == null)
			return;
		
		String groupname = params.get("groupname").getAsString();
		String message = params.get("message").getAsString();
		
		//get email addresses
		HashMap<String,String> emails = new HashMap<String,String>();
		try {
			emails = dba().getGroupEmails(user.username, groupname);
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "send-group-message");
		if (emails == null)
		{
			msg.addProperty("status", "fail");
		}
		else
		{
			for (String u : emails.keySet())
			{
				send_mail.enqueueMessage("Update for "+groupname,
						emails.get(u), message);
			}
			msg.addProperty("status", "success");
		}
		user.sendMessage(msg.toString());
	}

	/**
	 * Changes the calling user's password, notifies the user of the change if successful.
	 * The user may continue to authenticate via the session resume feature if their cookie
	 * data is preserved.
	 * @param s - the websocket session
	 * @param params - old (old password), password1 (new password, first repetition),
	 * password2 (new password, second repetition)
	 */
	public void changeUserPassword(Session s, JsonObject params) 
	{
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		if (params.get("old") == null || params.get("password1") == null || params.get("password2") == null)
			return;
		
		String old = params.get("old").getAsString();
		String password1 = params.get("password1").getAsString();
		String password2 = params.get("password2").getAsString();
		
		if (!password1.equals(password2) || old.equals(password1))
		{
			System.out.println("Passwords do not match, or did not change.");
			return;
		}

		old = this.sha512(old);
		password1 = this.sha512(password1);
		
		params.addProperty("rtype", "change-password");
		params.remove("old");
		params.remove("password1");
		params.remove("password2");
		
		try {
			if (dba().changeUserPassword(user.username, old, password1))
			{
				params.addProperty("status", "success");
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
		user.sendMessage(params.toString());
	}
	
	/**
	 * Generates a new non-specific authentication token by hashing the output of 
	 * a randomly generated number sequence concatenated with the user's username.
	 * @param user - uses the username to ensure the token will be unique to them
	 * @return the token as a 128 character hex string
	 */
	private String generateToken(String username)
	{
		//generate a token
		String pre = "";
		for (int i=0; i<128; i++)
		{
			pre += ""+sr.nextLong();
		}
		pre += username;
		String token = this.sha512(pre);
		
		return token;
	}

	/**
	 * Invokes a token request to change the email address of the calling user. Generates
	 * a new token and emails it to the specified email address. The user must retrieve the token
	 * to authorise the change.
	 * @param s - the websocket session
	 * @param params - email (the new email address)
	 */
	public void changeUserEmail(Session s, JsonObject params) 
	{
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		if (params.get("email") == null)
			return;
		
		String email = params.get("email").getAsString();
		params.addProperty("rtype", "change-email");
		
		String token = this.generateToken(user.username).substring(0,32);
		
		send_mail.enqueueMessage("Email Change Request", email, 
				"You have requested to change your primary email address to this one.<br>"
				+ "Follow this link to authorise this change:<br><br>"
				+ DOMAIN+"?token="+token+"<br><br>"
				+ "Regards,<br>"
				+ "PBrowse Team.");
		
		JsonObject tokparams = new JsonObject();
		tokparams.addProperty("email", email);
		
		try {
			if (dba().setToken(user.username, token, DBA.FUNCTION_CHANGEEMAIL, tokparams.toString()))
			{
				params.addProperty("status", "success");
			}
		} catch (SQLException e1) {
			e1.printStackTrace();
		}
		
		user.sendMessage(params.toString());
	}
	
	/**
	 * Invokes a token request to reset the password of the calling user. Generates
	 * a new token and emails it to the specified email address. The user must retrieve the token
	 * to authorise the change.
	 * @param s - the websocket session
	 * @param params - email (the new email address)
	 */
	public void resetUserPassword(Session s, JsonObject params) 
	{
		if (params.get("usermail") == null)
			return;
		
		params.addProperty("rtype", "reset-password");
		params.addProperty("status", "fail");
		
		String usermail = params.get("usermail").getAsString();
		HashMap<String, String> r;
		try {
			r = dba().lookupUser(usermail);
			if (r != null)
			{
				String email = r.get("email");
				String username = r.get("username");
				
				String token = this.generateToken(username).substring(0,32);
				
				send_mail.enqueueMessage("Password Reset Request", email, 
						"A password reset has been requested for this account. If you have not made this "
								+"request, simply ignore this message.<br>"
								+ "Follow this link to authorise this change:<br><br>"
								+ DOMAIN+"?token="+token+"<br><br>"
								+ "Regards,<br>"
								+ "PBrowse Team.");
				
				JsonObject tokparams = new JsonObject();
				tokparams.addProperty("newpassword", this.generateSID(12));
				
				if (dba().setToken(username, token, DBA.FUNCTION_RESETPASSWORD, tokparams.toString()))
				{
					params.addProperty("status", "success");
				}
			}
			
		} catch (SQLException e2) {
			//noop
		}
		
		try
		{
			s.getRemote().sendStringByFuture(params.toString());
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	/**
	 * Changes the nickname of the calling user and notifies them of its success
	 * @param s - the websocket session
	 * @param params - nickname
	 */
	public void changeUserNickname(Session s, JsonObject params) 
	{
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		if (params.get("nickname") == null)
			return;
		
		String nickname = params.get("nickname").getAsString();
		
		params.addProperty("rtype", "change-nickname");
		params.addProperty("username", user.username);

		try {
			if (dba().changeUserNickname(user.username, nickname))
			{
				params.addProperty("status", "success");
				user.nick = nickname;
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
		user.sendMessage(params.toString());
	}

	private HashMap<IUser,String> t_usergroup = new HashMap<IUser,String>();
	private HashMap<String,ArrayList<IUser>> t_groupusers = new HashMap<String,ArrayList<IUser>>();
	
	/**
	 * Message handler for all TEST functions i.e. usability test synchronisation
	 * @param s - the websocket session
	 * @param type 
	 * @param params - variable
	 */
	public void testSync(Session s, String type, JsonObject params) 
	{
		IUser user = ws_user.get(s);
		ArrayList<IUser> ul = null;
		
		switch (type)
		{
		case "TEST-set-group":
			String group = params.get("groupname").getAsString();
			
			if (t_usergroup.containsKey(user))
			{
				return;
			}
			t_usergroup.put(user, group);
			
			if (!t_groupusers.containsKey(group))
			{
				ul = new ArrayList<IUser>();
			}
			else
			{
				ul = t_groupusers.get(group);
			}
			ul.add(user);
			t_groupusers.put(group, ul);
			break;
		case "TEST-begin-question":
			params.addProperty("rtype",type);
			//get users in this caller's group to forward message to
			if (t_usergroup.containsKey(user))
			{
				ul = t_groupusers.get(t_usergroup.get(user));
				for (IUser u : ul)
				{
					if (u.equals(user)) continue;
					
					u.sendMessage(params.toString());
				}
			}
			break;
		}
		
	}
	
	/**
	 * Allows a current session leader to elect a new leader, without leaving the session. Additionally,
	 * the elected leader will retain this position if they leave and rejoin the session
	 * @param s - the WebSocket session
	 * @param params - leader (the new leader)
	 */
	public void setSessionLeader(Session s, JsonObject params)
	{
		params.addProperty("rtype","session-nominate-leader");
		
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		if (params.get("leader") == null)
			return;
		
		String leader = params.get("leader").getAsString();
		IUser target = username_user.get(leader);
		if (target == null)
		{
			params.addProperty("status","fail");
			params.addProperty("error","The nominated user is not available.");
			user.sendMessage(params.toString());
			return;
		}
		
		PSession ls = getUserSession(user);
		PSession us = getUserSession(target);
		if (!ls.equals(us) || ls == null)
		{
			params.addProperty("status","fail");
			params.addProperty("error","The nominated user is not part of this session.");
			user.sendMessage(params.toString());
			return;
		}
		
		if (!user.isLeader || target.isLeader)
		{
			params.addProperty("status","fail");
			params.addProperty("error","Cannot assign leadership.");
			user.sendMessage(params.toString());
			return;
		}
		
		params.addProperty("status","success");
		params.addProperty("caller",user.username);
		
		//do update
		user.isLeader = false; target.isLeader = true;
		//sets creator to new leader
		ls.creator = target;
		
		//send out updated statuses
		params.add("user-status", ls.getUsersStatus());
		ls.updateAllUsers(params.toString());
	}
	
	/**
	 * Adds or removes a user from the session blacklist - does not perform checks on the
	 * validity of the provided username, and will always succeed
	 * @param s - the WebSocket session
	 * @param params - user (the user to blacklist), option ("allow" or "deny")
	 */
	public void setSessionBlacklist(Session s, JsonObject params)
	{
		params.addProperty("rtype","session-blacklist-update");
		params.addProperty("status","success");
		
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		if (params.get("user") == null || params.get("option") == null)
			return;
		
		String b_user = params.get("user").getAsString();
		String option = params.get("option").getAsString();
		
		PSession ps = this.getUserSession(user);
		if (ps == null)
			return;
		
		if (!user.isLeader)
			return;
		
		//do blacklisting
		if (!ps.blacklistUser(b_user, option))
		{
			params.addProperty("status","fail");
			user.sendMessage(params.toString());
			return;
		}
		
		//successful blacklisting
		ps.updateAllUsers(params.toString());
		
		//if blacklisted user is in our session, automatically kick him
		if (option.equals("deny"))
			this.doSessionKickUser(s, params);
	}
	
	/**
	 * Invites a user to the current session, they must not already be in a session of
	 * their own, and the inviter must be leader. The invited user can always refuse the
	 * invitation.
	 * @param s - the WebSocket session
	 * @param params - user (the user to invite)
	 */
	public void doSessionInviteUser(Session s, JsonObject params)
	{
		params.addProperty("rtype","session-invite-user");
		params.addProperty("status","success");
		
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		if (params.get("user") == null)
			return;
		
		PSession ps = this.getUserSession(user);
		if (ps == null)
			return;
		
		
		String new_user = params.get("user").getAsString();
		IUser n_user = username_user.get(new_user);
		if (n_user == null)
		{
			params.addProperty("status","fail");
			params.addProperty("error","User not available.");
			user.sendMessage(params.toString());
			return;
		}
		
		//user already in another session
		PSession ns = this.getUserSession(n_user);
		if (ns != null)
		{
			params.addProperty("status","fail");
			params.addProperty("error","User cannot be invited.");
			user.sendMessage(params.toString());
			return;
		}
		
		params.addProperty("uuid",ps.getUUID());
		//send code if required
		String code = ps.getOption("code");
		if (code != null)
		{
			params.addProperty("code", code);
		}
		
		n_user.sendMessage(params.toString());
		user.sendMessage(params.toString());
	}
	
	/**
	 * Removes the specified user from the session, but does not prevent them from returning,
	 * must be called by session leader
	 * @param s - the WebSocket session
	 * @param params - user (the user to be removed)
	 */
	public void doSessionKickUser(Session s, JsonObject params)
	{
		params.addProperty("rtype","session-kick-user");
		params.addProperty("status","success");
		
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		if (params.get("user") == null)
			return;
		
		String new_user = params.get("user").getAsString();
		IUser target = username_user.get(new_user);
		
		//get user sessions - they should be the same
		PSession ps = this.getUserSession(user);
		PSession ns = this.getUserSession(target);
		if (ps == null || ns == null)
			return;
		
		if (!user.isLeader || target.isLeader)
		{
			return;
		}

		//cannot kick user in different session
		if (!ps.equals(ns))
		{
			return;
		}
		
		user.sendMessage(params.toString());
		target.sendMessage(params.toString());
	}
	
	/**
	 * Register a remote data source for usage in PBrowse - stores metadata only
	 * @param s
	 * @param params
	 */
	public void registerRemoteSource(Session s, JsonObject params)
	{
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		if (params.get("dataurl") == null || params.get("data_format") == null || 
				params.get("feature_seq") == null || params.get("studyid") == null || 
				params.get("trackname") == null || params.get("genometag") == null || 
				params.get("description") == null || params.get("ispublic") == null)
		{
			//missing parameters, noop
			return;
		}
		
		String dataurl = params.get("dataurl").getAsString();
		String data_format = params.get("data_format").getAsString();
		String feature_seq = params.get("feature_seq").getAsString();
		String studyid = params.get("studyid").getAsString();
		String trackname = params.get("trackname").getAsString();
		String genometag = params.get("genometag").getAsString();
		String description = params.get("description").getAsString();
		int ispublic = params.get("ispublic").getAsInt();
		
		//description includes genometag
		description = genometag+(description.length()>0?": "+description:"");
		
		//format will be used as a meta-information string since
		//it doesn't serve a crucial purpose anyway
		JsonObject meta = new JsonObject();
		meta.addProperty("studyid", studyid);
		meta.addProperty("data_format", data_format);
		meta.addProperty("feature_seq", feature_seq);
		
		String [] p = dataurl.split("/");
		String filename = p[p.length-1];
		meta.addProperty("filename", filename);
		
		DataDesc dd = null;
		int fid = -1;
		try {
			fid = dba().registerDataFile(user.username, //user who registers the file 'owns' it
					description, trackname, //standard
					dataurl, meta.toString(), //path becomes URL to remote file
					ispublic, null, //don't need index
					false, true //is a remotely located file
					);
			if (fid > 0)
			{
				//fetch as dd so we can return the normal parameters
				dd = dba().getDataFile(fid);
			}
		} catch (SQLException e) {
			e.printStackTrace();
			return;
		}
		
		if (dd != null)
		{
			JsonObject msg = dd.getJson();
			msg.addProperty("rtype","new-session-file");
			msg.addProperty("status","success");
			user.sendMessage(msg.toString());
		}
	}

	/**
	 * Stores the user provided configuration set for recall at a later date
	 * @param session - the WebSocket Session
	 * @param params - the config list
	 */
	public void updateSavedConfig(Session s, JsonObject params) 
	{
		IUser user = ws_user.get(s);
		if (!user.loggedIn)
		{
			return;
		}
		
		if (params.get("config") == null)
		{
			return;
		}
		String config = params.get("config").getAsString();
		
		boolean success = false;
		try {
			success = dba().updateStoredConfig(user.username, config);
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
		JsonObject msg = new JsonObject();
		msg.addProperty("rtype", "update-saved-config");
		msg.addProperty("status", (success?"success":"fail"));
		user.sendMessage(msg.toString());
	}
}
