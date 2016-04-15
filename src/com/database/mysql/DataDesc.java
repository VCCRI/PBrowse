package com.database.mysql;

import com.google.gson.JsonObject;

/**
 * Data Descriptor object, used in recalling file metadata from the database
 * @author root
 *
 */
public class DataDesc 
{
	public long id;
	public int access_required;
	public String owner;
	public String description;
	public String path;
	public String format;
	public String trackname;
	public int ispublic;
	public Long indexedby;
	public boolean isindex;
	public JsonObject groups; //groupname : {required_access, user_counts}
	public boolean remote;
	
	public DataDesc(
			long id, String owner, int access_required,
			String description, String path, String trackname, 
			String format, int ispublic, long indexedby, boolean isindex,
			JsonObject groups, boolean remote
			)
	{
		this.id = id;
		this.owner = owner;
		this.access_required = access_required;
		this.description = description;
		this.path = path;
		this.format = format;
		this.trackname = trackname;
		this.ispublic = ispublic;
		this.indexedby = indexedby;
		this.isindex = isindex;
		this.groups = groups;
		this.remote = remote;
	}
	
	/**
	 * Gets all of the DD fields as a Json object for easy transmission to the client
	 * @return
	 */
	public JsonObject getJson()
	{
		JsonObject me = new JsonObject();
		me.addProperty("id", id);
		me.addProperty("owner", owner);
		me.addProperty("description", description);
		me.addProperty("path", path);
		me.addProperty("format", format);
		me.addProperty("trackname", trackname);
		me.addProperty("ispublic", ispublic);
		me.addProperty("indexedby", indexedby);
		me.addProperty("isindex", isindex);
		me.addProperty("access_required", access_required);
		me.addProperty("remote", remote);
		return me;
	}
}
