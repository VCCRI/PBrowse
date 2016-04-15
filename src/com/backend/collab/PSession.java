package com.backend.collab;

import java.util.HashSet;
import java.util.Set;
import java.util.HashMap;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonPrimitive;

/**
 * Session manager manages instances of this object, which contains grouped
 * user instances
 *
 * Abstract representation of a collaborative session
 *
 * PSession belongs to SessionManager
 * @author root
 *
 */

public class PSession 
{
	//set of users in the session
	private final HashSet<IUser> users = new HashSet<IUser>();
	private String uuid = null;

	//mapping of option parameters
	private final HashMap<String,String> options = new HashMap<String,String>();
	
	//the current set of tiers viewed by the session users
	private JsonArray tiers = null;
	
	//persistent set of joining user names, allows for sharing of private files in the session
	private final HashSet<String> owners = new HashSet<String>();
	
	//storage of chat history for the session
	public final JsonArray chatHistory = new JsonArray();
	
	//list of users blocked from joining the session username -> null
	private final HashSet<String> blacklist = new HashSet<String>();
	
	public IUser creator = null;
	public String genome = null;
	public String custom_genome = null;
	
	//return jsonobject of options
	public JsonObject getOptions()
	{
		JsonObject opts = new JsonObject();
		for (String opt : options.keySet())
		{
			opts.addProperty(opt, options.get(opt));
		}
		return opts;
	}
	//setter for options
	public void setOption(String opt, String val)
	{
		options.put(opt,val);
	}
	
	public JsonArray getTiers()
	{
		return this.tiers;
	}
	public void setTiers(JsonArray t)
	{
		this.tiers = t;
	}
	
	/**
	 * Gets the value of a particular option
	 * @param opt
	 * @return
	 */
	public String getOption(String opt)
	{
		if (options.containsKey(opt))
			return options.get(opt);
		
		else return null;
	}

	/**
	 * create a new session with provided first user as leader
	 * @param user
	 */
	public PSession(IUser user)
	{
		//leader
		user.isLeader = true;
		creator = user;
		users.add(user);
		owners.add(user.username);

		//default options
		options.put("allowdivergence","true");
		options.put("followersync","true");
		options.put("private","true");
	}
	public void setSID(String uuid)
	{
		this.uuid = uuid;
	}
	public final String getSID()
	{
		return this.uuid;
	}
	
	/**
	 * Adds a new user to the session, returns false if the user is blacklisted
	 * @param user
	 * @return
	 */
	public boolean addUser(IUser user)
	{
		if (blacklist.contains(user.username))
			return false;
		
		if (user.username == creator.username)
		{
			//demote current leader
			for (IUser us : users)
			{
				if (us.isLeader)
				{
					us.isLeader = false;
					break;
				}
			}
			user.isLeader = true;
		}
		else
		{
			user.isLeader = false;
		}
		users.add(user);
		//register owner
		owners.add(user.username);
		return true;
	}
	
	public boolean hasAccess(String owner)
	{
		return owners.contains(owner);
	}
	
	public HashSet<IUser> userlist()
	{
		return users;
	}
	
	/**
	 * Removes the specified user from the session
	 * @param user - user to remove
	 * @return number of users remaining in the session
	 */
	public int removeUser(IUser user)
	{
		if (user.isLeader)
		{
			//promote another user to leader
			for (IUser u : users)
			{
				if (u.equals(user))
				{
					continue;
				}
				else
				{
					//promote the next user to leader status
					u.isLeader = true;
					break;
				}
			}
			
		}
		
		users.remove(user);

		//if there are no other users to promote, remove session entirely
		return users.size();
	}
	
	public void closeSession()
	{
		users.clear();
	}
	
	public String getUUID()
	{
		return uuid;
	}
	
	/**
	 * Gets an array of all the users statuses
	 * @return
	 */
	public JsonArray getUsersStatus()
	{
		JsonArray status = new JsonArray();
		
		for (IUser user : users)
		{
			status.add(user.getStatus());
		}
		
		return status;
	}
	
	/**
	 * Get the session blacklist as a JSon array to return to new users joining the session
	 * @return
	 */
	public JsonArray getBlacklist()
	{
		JsonArray bl = new JsonArray();
		
		for (String str : this.blacklist)
		{
			bl.add(new JsonPrimitive(str));
		}
		
		return bl;
	}
	
	public int numUsers()
	{
		return users.size();
	}
	
	public Set<String> getOwners()
	{
		return owners;
	}
	
	/**
	 * Sends the given message string to all users in the session
	 * For updating changes to users
	 * @param msg
	 */
	public void updateAllUsers(String msg)
	{
		for (IUser u : users)
		{
			if (u == null)
				continue;
			
			u.sendMessage(msg);
		}
	}
	
	public boolean equals(PSession p)
	{
		if (this.uuid == p.uuid)
		{
			return true;
		}
		return false;
	}
	
	/**
	 * adds or removes a user from the session blacklist
	 * @param user - username of the user
	 * @param opt - allow, deny
	 * @return success
	 */
	public boolean blacklistUser(String user, String opt)
	{
		switch (opt)
		{
		case "allow":
			blacklist.remove(user);
			return true;
		case "deny":
			blacklist.add(user);
			return true;
		}
		return false;
	}
	
}