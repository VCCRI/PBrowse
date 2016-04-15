package com.database.mysql;

import com.google.gson.JsonObject;

/**
 * Persistent Comment/bookmark object
 * @author peter
 *
 */

public class Comment 
{
	public long id;
	public long trackref;
	public String ctext;
	public String chr;
	public long startpos;
	public long endpos;
	public String creator;
	public int ispublic;
	
	public Comment(long id, long trackref, String ctext, String chr, long startpos, long endpos, String creator, int ispublic)
	{
		this.id = id;
		this.trackref = trackref;
		this.ctext = ctext;
		this.chr = chr;
		this.startpos = startpos;
		this.endpos = endpos;
		this.creator = creator;
		this.ispublic = ispublic;
	}
	
	/**
	 * Gets the field set as a Json Object for direct transmission to the client
	 * @return
	 */
	public JsonObject getJson()
	{
		JsonObject me = new JsonObject();
		me.addProperty("id", id);
		me.addProperty("trackref", trackref);
		me.addProperty("ctext", ctext);
		me.addProperty("chr", chr);
		me.addProperty("startpos", startpos);
		me.addProperty("endpos", endpos);
		me.addProperty("creator", creator);
		me.addProperty("ispublic", ispublic);
		return me;
	}
}
