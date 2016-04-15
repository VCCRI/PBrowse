package com.backend.collab;

import org.eclipse.jetty.websocket.api.Session;

import com.google.gson.JsonObject;

/**
 * The user data which is maintained by the server for all connected
 * users. Updates to this object will trigger broadcasts to all users
 * within the same interactive session
 *
 * IUser belong to PSession, may also be managed by SessionManager
 * @author root
 *
 */

public class IUser 
{
	//the websocket session used to communicate with this user
	private Session session;
	
	/**
	 * this is the information which will be passed around and updated
	 */
	public String chr = "chr1";
	public long start_win = 0;
	public long end_win = 0;

	public String username = null;
	public String nick = null;
	private String uauth = null;
	public boolean loggedIn = false;
	
	//groupname : {access_level : x, counts : y }
	public JsonObject membership = null;

	//determines whether a user controls the session
	public boolean isLeader;
	
	/**
	 * Clears all user information from this object
	 */
	public void clearData()
	{
		username = null;
		nick = null;
		uauth = null;
		loggedIn = false;
		membership = null;
	}
	
	public IUser(Session s)
	{
		session = s;
		isLeader = false;
	}
	
	public IUser(Session s, boolean isleader)
	{
		session = s;
		isLeader = isleader;
	}
	
	/**
	 * Changes the uauth token
	 * @param a
	 */
	public void setUAuth(String a)
	{
		if (uauth == null)
		{
			uauth = a;
		}
	}
	
	public String getUAuth()
	{
		return uauth;
	}
	
	/**
	 * packages status into json object for transmission
	 * @return
	 */
	public JsonObject getStatus()
	{
		JsonObject s = new JsonObject();
		s.addProperty("chr", chr);
		s.addProperty("start_win", start_win);
		s.addProperty("end_win", end_win);
		s.addProperty("username", username);
		s.addProperty("nickname", nick);
		s.addProperty("isLeader", isLeader);
//		s.add("groups", membership);
		return s;
	}
	
	/**
	 * Updates the server's knowledge of the user's position
	 * @param params - positional information
	 */
	public void updateUserPos(JsonObject params)
	{
		chr = params.get("chr").getAsString();
		start_win = params.get("start_win").getAsLong();
		end_win = params.get("end_win").getAsLong();
	}
	
	/**
	 * Sends a message to this user
	 * @param message - the JSON encoded message
	 */
	public void sendMessage(String message)
	{
		if (session == null) return;
		
		try
		{
			session.getRemote().sendStringByFuture(message);
		}
		catch (Exception e)
		{
			//the socket was closed
			session = null;
		}
	}

	public void setSession(Session s) 
	{
		this.session = s;
	}
	
	public Session getSession()
	{
		return this.session;
	}
}