package com.database.mysql;

import java.io.File;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Types;
import java.util.ArrayList;
import java.util.HashMap;

import net.sf.ehcache.Element;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.mysql.jdbc.jdbc2.optional.MysqlConnectionPoolDataSource;

import static com.backend.collab.CacheManager.cm;

/**
 * MySQL adapter for DB calls
 * @author peter
 *
 */
public class DBA 
{
	private static DBA dba = null;
	private MysqlConnectionPoolDataSource dataSource = null;
	
	/**
	 * Permissions for group access
	 */
	public static final int PERMISSION_OWNERSHIP 	= 0b1000;
	public static final int PERMISSION_MANAGEUSERS 	= 0b0100;
	public static final int PERMISSION_MANAGEFILES	= 0b0010;
	public static final int PERMISSION_READ 		= 0b0001;
	
	/**
	 * Token functions
	 */
	public static final int FUNCTION_RESETPASSWORD   = 0b100;
	public static final int FUNCTION_REGISTERACCOUNT = 0b010;
	public static final int FUNCTION_CHANGEEMAIL	 = 0b001;
	
	/**
	 * Privacy Constants
	 */
	public static final int PRIVACY_ONLYME 	= 0b001;
	public static final int PRIVACY_PRIVATE	= 0b010;
	public static final int PRIVACY_PUBLIC	= 0b100;
	
	/**
	 * The DBA is a singleton instance
	 */
	private DBA()
	{
		dataSource = new MysqlConnectionPoolDataSource();
		dataSource.setUser("pbrowse");
		dataSource.setPassword("pbrowse");
		dataSource.setServerName("localhost");
	}
	
	public static DBA dba()
	{
		if (dba == null)
		{
			dba = new DBA();
		}
		return dba;
	}
	
	/**
	 * Permission verifier, if the provided acl has permission p, return true
	 * @param acl the user's privilege
	 * @param p the privilege to check against
	 * @return true if acl has privilege p
	 */
	public static boolean P(int acl, int p)
	{
		if ((acl&p) == p) return true;
		return false;
	}
	
	/**
	 * Verify the creds of the provided user, allowing them to log in to the system
	 * 
	 * @param username - username of logging in user
	 * @param password - hashed password string
	 * @return user data or null for bad login
	 * @throws SQLException 
	 */
	public HashMap<String,String> doLogin(String username, String password) throws SQLException
	{
		if (username == null || password == null)
		{
			System.err.println("Null values provided in login.");
			return null;
		}
		
		PreparedStatement s = null;
		ResultSet r = null;
		HashMap<String,String> data = null;
		Connection con = null;
		try {
			con = dataSource.getConnection();
			
			s = con.prepareStatement("SELECT username, nick, activated "
					+ "FROM pbrowsedb.users "
					+ "WHERE username = ? AND password = ?");
			s.setString(1, username);
			s.setString(2, password);
			
			r = s.executeQuery();
			if (r.next())
			{
				if (r.getBoolean(3) == true)
				{
					data = new HashMap<String,String>();
					data.put("username", r.getString(1));
					data.put("nick", r.getString(2));
				}
			}
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
		}
		finally 
		{
			if (con != null) con.close();
			if (r != null) r.close();
			if (s != null) s.close();
		}
		return data;
	}
	
	public String getUserSavedConfig(String username) throws SQLException
	{
		if (username == null)
		{
			return null;
		}
		
		PreparedStatement s = null;
		ResultSet r = null;
		Connection con = null;
		
		String result = null;
		try {
			con = dataSource.getConnection();
			
			s = con.prepareStatement("SELECT saved_config "
					+ "FROM pbrowsedb.users "
					+ "WHERE username = ?;");
			s.setString(1, username);
			
			r = s.executeQuery();
			if (r.next())
			{
				result = r.getString(1);
			}
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
			result = null;
		}
		finally 
		{
			if (con != null) con.close();
			if (r != null) r.close();
			if (s != null) s.close();
		}
		return result;
	}
	
	/**
	 * Registers a new user account for use with PBrowse
	 * @param username - the unique username
	 * @param nick - the pseudonym of the user
	 * @param email - the user's email
	 * @param password - the user's password
	 * @return success if the username is not already taken
	 * @throws SQLException
	 */
	public boolean doRegister(String username, String nick, String email, String password) throws SQLException
	{
		if (username == null || password == null || nick == null || email == null)
		{
			System.err.println("Null values provided in registration.");
			return false;
		}
		
		PreparedStatement s = null;
		boolean result = false;
		Connection con = null;
		try {
			con = dataSource.getConnection();
			
			s = con.prepareStatement("INSERT INTO pbrowsedb.users (username,email,password,nick)"
					+ "VALUES(?,?,?,?);");
			s.setString(1, username);
			s.setString(2, email);
			s.setString(3, password);
			s.setString(4, nick);
			
			s.executeUpdate();
			s.close();
			con.close();
			result = true;
		} 
		catch (SQLException e) 
		{
			System.out.println("Username is already taken.");
			result = false;
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
		}
		return result;
	}
	
	/**
	 * Gets list of all data files belonging to the user
	 * @param username - the calling user
	 * @return array of files as Data Descriptors
	 * @throws SQLException 
	 */
	@SuppressWarnings("unchecked")
	public ArrayList<DataDesc> getDataFiles(String username) throws SQLException
	{
		if (username == null)
		{
			System.err.println("No user provided.");
			return null;
		} 
		
		ArrayList<DataDesc> data = null;
		
		//attempt to recall file list from cache to minimize DB calls
        if (cm().userfilesCache.isKeyInCache(username))
        {
        	Element e = cm().userfilesCache.get(username);
        	if (e != null)
        	{
	        	if (!e.isExpired())
	        	{
        			data = (ArrayList<DataDesc>) e.getObjectValue();
	        	}
        	}
        }

        //the cache entry expired or does not exist
        if (data == null)
        {
        	PreparedStatement s = null;
    		ResultSet r = null;
    		Connection con = null;
    		try {
    			con = dataSource.getConnection();
    			
    			s = con.prepareStatement("SELECT id, owner, description, path, trackname, format, ispublic, indexedby, remote "
    					+ "FROM pbrowsedb.data "
    					+ "WHERE owner = ? AND isindex = 0");
    			s.setString(1, username);
    			
    			data = new ArrayList<DataDesc>();
    			r = s.executeQuery();
    			while (r.next())
    			{
    				//filter out onlyme private files if the caller is not the owner
    				if (r.getInt(7) == DBA.PRIVACY_ONLYME)
    				{
    					if (!r.getString(2).equals(username))
    						continue;
    				}
    				
    				//not retrieving index files
    				data.add( new DataDesc(r.getLong(1),
    						r.getString(2),
    						0, //access_required, irrelevant here
    						r.getString(3),
    						r.getString(4),
    						r.getString(5),
    						r.getString(6),
    						r.getInt(7),
    						r.getLong(8), 
    						false, //is file an index, only return non-index entries
    						null, //group list
    						r.getBoolean(9)) //is a remote file
    				);
    			}
    			//add entry to cache
    			cm().userfilesCache.put(new Element(username,data));
    		} 
    		catch (SQLException e) 
    		{
    			e.printStackTrace();
    			data = null;
    		}
    		finally 
    		{
    			if (con != null) con.close();
    			if (r != null) r.close();
    			if (s != null) s.close();
    		}
        }
		
		return data;
	}
	
	/**
	 * Gets data file with the provided id - returns a cached entry if the data was not modified
	 * since the last request for it
	 * @param username - the calling user
	 * @return dd of stored file
	 * @throws SQLException 
	 */
	public DataDesc getDataFile(long fileid) throws SQLException
	{
		PreparedStatement s = null;
		ResultSet r = null;
		DataDesc data = null;
		Connection con = null;
		
		//attempt to recall DD from cache to minimize DB calls
        if (cm().ddCache.isKeyInCache(fileid))
        {
        	Element e = cm().ddCache.get(fileid);
        	if (e != null)
        	{
	        	if (!e.isExpired())
	        	{
	        		data = (DataDesc) e.getObjectValue();
	        	}
        	}
        }

        //the cache entry expired or does not exist
        if (data == null)
        {
        	try {
        		con = dataSource.getConnection();
        		
        		System.out.println("requested: "+fileid);
        		s = con.prepareStatement("SELECT id, owner, description, path, trackname, format, ispublic, indexedby, isindex, remote "
        				+ "FROM pbrowsedb.data "
        				+ "WHERE id = ?");
        		s.setLong(1, fileid);
        		
        		//get access groups
        		JsonObject groups = this.getFileGroups(fileid);
        		
        		r = s.executeQuery();
        		if (r.next())
        		{
        			data = new DataDesc(
        					fileid,r.getString(2),
        					-1,r.getString(3),r.getString(4),
        					r.getString(5),r.getString(6),
        					r.getInt(7),r.getLong(8),
        					r.getBoolean(9),groups, r.getBoolean(10)
        					);
        		}
        	} 
        	catch (SQLException e) 
        	{
        		e.printStackTrace();
        		data = null;
        	}
        	finally 
        	{
        		if (con != null) con.close();
        		if (r != null) r.close();
        		if (s != null) s.close();
        	}
        	//otherwise cache the new entry
        	cm().ddCache.put(new Element(fileid,data));
        }
		return data;
	}
	
	/**
	 * Stores the Data Descriptor of the new file
	 * @param username - owner of the file
	 * @param desc - long description of the data
	 * @param trackname - identifying name when loaded in the browser
	 * @param path - the location of the file on the filesystem
	 * @param format - the type of data
	 * @return id of the registered data file
	 * @throws SQLException 
	 */
	public int registerDataFile(
			String username, String desc, 
			String trackname, String path, 
			String format, int ispublic, 
			Integer indexid, boolean isindex, boolean remote) throws SQLException
	{
		if (username == null)
		{
			System.err.println("No user provided.");
			return -1;
		}
		
		PreparedStatement s = null;
		PreparedStatement s2 = null;
		ResultSet r = null;
		int result = 0;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			
			s = con.prepareStatement("INSERT INTO pbrowsedb.data "
					+ "(owner,description,trackname,path,format,ispublic,indexedby,isindex,remote) "
					+ "VALUES(?,?,?,?,?,?,?,?,?)");
			s.setString(1, username);
			s.setString(2, desc);
			s.setString(3, trackname);
			s.setString(4, path);
			s.setString(5, format);
			s.setInt(6, ispublic);
			
			if (indexid == null) s.setNull(7, Types.NULL);
			else s.setInt(7, indexid);
			
			s.setBoolean(8, isindex);
			s.setBoolean(9, remote);
			
			s.executeUpdate();
			
			s2 = con.prepareStatement("SELECT LAST_INSERT_ID();");
			r = s2.executeQuery();
			if (r.next())
			{
				result = r.getInt(1);
			}
			
			//invalidate user files list only
			cm().userfilesCache.remove(username);
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
			result = -1;
		}
		finally 
		{
			if (con != null) con.close();
			if (r != null) r.close();
			if (s != null) s.close();
			if (s2 != null) s2.close();
		}
		return result;
	}
	
	/**
	 * Removes the entry for the data descriptor with given id
	 * @param username - the file owner
	 * @param id - the id of the file
	 * @return success if the calling user owns the file
	 * @throws SQLException
	 */
	public boolean deleteDataFile(String username, long id) throws SQLException
	{
		if (username == null)
		{
			System.err.println("No user provided.");
			return false;
		}
		
		PreparedStatement s = null;
		PreparedStatement s2 = null;
		ResultSet r = null;
		boolean success = true;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			
			s = con.prepareStatement("SELECT id, path, remote FROM pbrowsedb.data WHERE (id = "
					+ "(SELECT indexedby FROM pbrowsedb.data WHERE id = ?) "
					+ "OR id = ?) AND owner = ?;");
			s.setLong(1, id);
			s.setLong(2, id);
			s.setString(3, username);
			
			long fileid = -1;
			String filepath = null;
			r = s.executeQuery();
			
			while (r.next())
			{
				fileid = r.getLong(1);
				filepath = r.getString(2);
				
				boolean remote = r.getBoolean(3);
				
				s2 = con.prepareStatement("DELETE FROM pbrowsedb.data "
						+ "WHERE owner = ? AND id = ?;");
				s2.setString(1, username);
				s2.setLong(2, fileid);
				
				if (!remote)
				{
					if (!(new File("/opt/pbrowse/sessions/"+username+filepath)).delete())
					{
						System.out.println("Failed to delete the file: "+"/opt/pbrowse/sessions/"+username+filepath);
						success = false;
						break;
					}
				}
				s2.executeUpdate();
				
				//invalidate
				cm().ddCache.remove(fileid);
				cm().userfilesCache.remove(username);
			}
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
			success = false;
		}
		finally 
		{
			if (con != null) con.close();
			if (r != null) r.close();
			if (s != null) s.close();
			if (s2 != null) s2.close();
		}
		
		return success;
	}
	
	/**
	 * Sets the state of a data file as being public or private
	 * @param username - the calling user
	 * @param id - file id
	 * @return success if the caller owns the file
	 * @throws SQLException
	 */
	public boolean setUserFilePublicStatus(String username, long id, int status) throws SQLException
	{
		if (username == null)
		{
			System.err.println("No user provided.");
			return false;
		}
		
		PreparedStatement s = null;
		PreparedStatement s2 = null;
		ResultSet r = null;
		boolean success = true;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			
			s = con.prepareStatement("SELECT indexedby "
					+ "FROM pbrowsedb.data "
					+ "WHERE owner = ? AND id = ?");
			s.setString(1, username);
			s.setLong(2, id);
			
			long indexid = 0;
			r = s.executeQuery();
			if (r.next())
			{
				indexid = r.getLong(1);
				s2 = con.prepareStatement("UPDATE pbrowsedb.data SET ispublic = ? "
						+ "WHERE owner = ? AND ( id = ? OR id = ? )");

				s2.setInt(1, status);
				s2.setString(2, username);
				s2.setLong(3, id);
				s2.setLong(4, indexid);
				
				s2.executeUpdate();
				
				//invalidate
				cm().ddCache.remove(id);
				cm().ddCache.remove(indexid);
				cm().userfilesCache.remove(username);
			}
			else success = false;
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
			success = false;
		}
		finally 
		{
			if (con != null) con.close();
			if (r != null) r.close();
			if (s != null) s.close();
			if (s2 != null) s2.close();
		}
		
		return success;
	}

	/**
	 * Gets a list of all the publically accessible files as Data Descriptors
	 * @return array of data descriptors
	 * @throws SQLException
	 */
	public ArrayList<DataDesc> getPublicFiles() throws SQLException
	{
		PreparedStatement s = null;
		ResultSet r = null;
		ArrayList<DataDesc> data = null;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			
			s = con.prepareStatement("SELECT id, owner, description, path, trackname, format, ispublic, indexedby, remote "
					+ "FROM pbrowsedb.data "
					+ "WHERE ispublic = ? AND isindex = 0");
			s.setInt(1, PRIVACY_PUBLIC);
			
			data = new ArrayList<DataDesc>();
			r = s.executeQuery();
			while (r.next())
			{
				data.add( new DataDesc(r.getLong(1),
						r.getString(2),
						0, //access required, irrelevant here
						r.getString(3),
						r.getString(4),
						r.getString(5),
						r.getString(6),
						r.getInt(7), 
						r.getLong(8), false, 
						null, r.getBoolean(9)) 
				);
			}
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
			data = null;
		}
		finally 
		{
			if (con != null) con.close();
			if (r != null) r.close();
			if (s != null) s.close();
		}
		
		return data;
	}
	
	/**
	 * Create a new comment based on parameters set in the provided Comment Object
	 * @param c - the Comment
	 * @return the id of the newly created comment
	 * @throws SQLException 
	 */
	public int newComment(Comment c) throws SQLException
	{
		PreparedStatement s = null;
		PreparedStatement s2 = null;
		ResultSet r = null;
		int result = 0;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			
			s = con.prepareStatement("INSERT INTO pbrowsedb.comments "
					+ "(trackref,ctext,chromosome,startpos,endpos,creator,ispublic) "
					+ "VALUES(?,?,?,?,?,?,?);");
			s.setLong(1, c.trackref);
			s.setString(2, c.ctext);
			s.setString(3, c.chr);
			s.setLong(4, c.startpos);
			s.setLong(5, c.endpos);
			s.setString(6, c.creator);
			s.setInt(7, c.ispublic);
			s.executeUpdate();
			
			s2 = con.prepareStatement("SELECT LAST_INSERT_ID();");
			r = s2.executeQuery();
			
			if (r.next())
			{
				result = r.getInt(1);
			}
			
			//invalidate cache
			cm().userCommentsCache.remove(c.trackref);
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
			result = 0;
		}
		finally 
		{
			if (con != null) con.close();
			if (r != null) r.close();
			if (s != null) s.close();
			if (s2 != null) s2.close();
		}
		
		return result;
	}
	
	/**
	 * Delete the comment with id
	 * @param cid - the id of the comment
	 * @param trackref 
	 * @param creator - the calling user
	 * @return success if the calling user owns the comment
	 * @throws SQLException
	 */
	public boolean deleteComment(long cid, long trackref, String creator) throws SQLException
	{
		PreparedStatement s = null;
		boolean success = true;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			
			s = con.prepareStatement("DELETE FROM pbrowsedb.comments WHERE id = ? AND trackref = ? AND creator = ?");
			s.setLong(1, cid);
			s.setLong(2, trackref);
			s.setString(3, creator);
			s.executeUpdate();
			s.close();
			con.close();
			
			//invalidate cache
			cm().userCommentsCache.remove(trackref);
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
			success = false;
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
		}
		return success;
	}
	
	/**
	 * Return all public comments for the given track, or all the comments made by the caller on the track
	 * @param <T> - any iterable set of username strings
	 * @param username - the calling user
	 * @param trackref - the id of the file the comment applies to
	 * @return array of comments
	 * @throws SQLException 
	 */
	@SuppressWarnings("unchecked")
	public <T extends Iterable<?>> ArrayList<Comment> getAllComments(T usernames, Long trackref, String caller) throws SQLException
	{
		ArrayList<Comment> comments = null;
		PreparedStatement s = null;
		ResultSet r = null;
		Connection con = null;
		
		//attempt to recall file list from cache to minimize DB calls
        if (cm().userCommentsCache.isKeyInCache(trackref))
        {
        	Element e = cm().userCommentsCache.get(trackref);
        	if (e != null)
        	{
	        	if (!e.isExpired())
	        	{
        			comments = (ArrayList<Comment>) e.getObjectValue();
	        	}
        	}
        }
		
        if (comments == null)
        {
        	try {
        		con = dataSource.getConnection();
        		
        		//fetch ALL comments on the track
        		s = con.prepareStatement("SELECT id, trackref, ctext, chromosome, startpos, endpos, creator, ispublic "
        				+ "FROM pbrowsedb.comments "
        				+ "WHERE trackref = ?;");
        		s.setLong(1, trackref);
        		
        		r = s.executeQuery();
        		comments = new ArrayList<Comment>();
        		while (r.next())
        		{
        			comments.add( new Comment(r.getLong(1),
        					r.getLong(2),
        					r.getString(3),
        					r.getString(4),
        					r.getLong(5),
        					r.getLong(6),
        					r.getString(7),
        					r.getInt(8)) 
        					);
        		}
        		//add entry to cache
    			cm().userCommentsCache.put(new Element(trackref,comments));
        	} 
        	catch (SQLException e) 
        	{
        		e.printStackTrace();
        		comments = null;
        	}
        	finally 
        	{
        		if (con != null) con.close();
        		if (r != null) r.close();
        		if (s != null) s.close();
        	}
        }
        
        //return comments which belong to users in the set of provided users
        ArrayList<Comment> filtered = new ArrayList<Comment>();
        for (Comment c : comments)
        {
        	if (c.ispublic == PRIVACY_PUBLIC)
        	{
        		filtered.add(c);
        	}
        	else if (c.ispublic == PRIVACY_ONLYME)
        	{
        		if (c.creator.equals(caller))
        		{
        			filtered.add(c);
        		}
        	}
        	else
        	{
        		for (Object uname : usernames)
        		{
        			if (c.creator.equals((String)uname))
        			{
        				filtered.add(c);
        			}
        		}
        	}
        }
        
		return filtered;
	}
	
	/**
	 * Gets all public comments or those made by the caller, which fit within the current view window
	 * @param username - calling user
	 * @param trackref - file id
	 * @param startpos - start genomic window
	 * @param endpos - ending genomic window
	 * @param chr - the chromosome
	 * @return array of comments
	 * @throws SQLException 
	 */
	public ArrayList<Comment> getCommentsInView(String username, long trackref, long startpos, long endpos, String chr) throws SQLException
	{
		ArrayList<Comment> comments = null;
		PreparedStatement s = null;
		ResultSet r = null;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			
			s = con.prepareStatement("SELECT id, trackref, ctext, chromosome, startpos, endpos, creator, ispublic "
					+ "FROM pbrowsedb.comments "
					+ "WHERE trackref = ? AND chromosome = ? AND (creator = ? OR ispublic = 1) AND "
					+ "( (startpos >= ? AND endpos <= ?) OR (startpos >= ? AND endpos <= ?) )");
			s.setLong(1, trackref);
			s.setString(2, chr);
			s.setString(3, username);
			s.setLong(4, startpos);
			s.setLong(5, startpos);
			s.setLong(6, endpos);
			s.setLong(7, endpos);
			
			r = s.executeQuery();
			comments = new ArrayList<Comment>();
			while (r.next())
			{
				comments.add( new Comment(r.getLong(1),
						r.getLong(2),
						r.getString(3),
						r.getString(4),
						r.getLong(5),
						r.getLong(6),
						r.getString(7),
						r.getInt(8) ) 
				);
			}
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
			comments = null;
		}
		finally 
		{
			if (con != null) con.close();
			if (r != null) r.close();
			if (s != null) s.close();
		}
		
		return comments;
	}
	
	/**
	 * Creates a new group with specified groupname. The calling user is
	 * automatically made owner of the group
	 * @param owner - the calling user
	 * @param groupname - the name of the new group
	 * @return success if the groupname is unique
	 * @throws SQLException
	 */
	public boolean createUserGroup(String owner, String groupname) throws SQLException
	{
		PreparedStatement s = null;
		Connection con = null;
		boolean success = false;
		
		try {
			con = dataSource.getConnection();
			
			s = con.prepareStatement("INSERT INTO pbrowsedb.groups "
					+ "(name, owner) "
					+ "VALUES(?,?);");
			s.setString(1, groupname);
			s.setString(2, owner);
			
			s.executeUpdate();
			
			//register that member as having owner level access and full permissions
			success = this.addUserToGroup(null, owner, groupname, 
					PERMISSION_OWNERSHIP
					|PERMISSION_MANAGEFILES
					|PERMISSION_MANAGEUSERS
					|PERMISSION_READ);
		} 
		catch (SQLException e) 
		{
			success = false;
			System.out.println("Group name already taken!");
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
		}
		
		return success;
	}
	
	/**
	 * Totally eradicates all group connections and file acls. ONLY the group creator can do this.
	 * @param caller - the calling user
	 * @param groupname - name of the group to delete
	 * @return success if the group was deleted
	 * @throws SQLException 
	 */
	public boolean deleteUserGroup(String caller, String groupname) throws SQLException
	{
		PreparedStatement s = null;
		boolean success = false;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			
			//verify ownership of caller
			int acl = this.getGroupAccessLevel(caller, groupname);
			
			if (P(acl,PERMISSION_OWNERSHIP))
			{
				s = con.prepareStatement("DELETE FROM pbrowsedb.groups "
						+ "WHERE name = ? AND owner = ?;");
				s.setString(1, groupname);
				s.setString(2, caller);
				
				s.executeUpdate();
				success = true;
				
				//invalidate
				cm().groupCache.remove(caller);
			}
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
			success = false;
		}
		finally
		{
			if (con != null) con.close();
			if (s != null) s.close();
		}
		return success;
	}
	
	/**
	 * Private function: returns the calling user's access level within the given group
	 * @param user - calling user
	 * @param group - the groupname
	 * @return access level
	 * @throws SQLException
	 */
	private int getGroupAccessLevel(String user, String group) throws SQLException
	{
		PreparedStatement s = null;
		ResultSet r = null;
		int acl = 0;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			s = con.prepareStatement("SELECT access_level "
					+ "FROM pbrowsedb.membership "
					+ "WHERE username = ? AND groupname = ?");
			s.setString(1, user);
			s.setString(2, group);
	
			r = s.executeQuery();
			if (r.next())
			{
				acl = r.getInt(1);
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
			acl = 0;
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
			if (r != null) r.close();
		}
		return acl;
	}
	
	/**
	 * Adds the specified 'user' to the group with the privileges of the 'caller' user
	 * @param caller - the user invoking the request
	 * @param user - the user to add to the group
	 * @param group - the name of the group
	 * @param access_level - the privilege level of the new user
	 * @return success if the user was added
	 * @throws SQLException
	 */
	public boolean addUserToGroup(String caller, String user, String group, int access_level) throws SQLException
	{
		PreparedStatement s = null;
		int acl = 0;
		boolean success = false;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			
			//the special case where this is the first user being added to the group
			//i.e. the creator, only the DBA can call this
			if (caller != null)
			{
				acl = getGroupAccessLevel(caller, group);
			}
			else
			{
				acl = PERMISSION_OWNERSHIP;
			}
			//only owners can set this permission
			if ( P(access_level,PERMISSION_MANAGEUSERS) && !P(acl,PERMISSION_OWNERSHIP) )
				success = false;
			//only pbrowse can set ownership
			else if ( P(access_level,PERMISSION_OWNERSHIP) && caller != null )
				success = false;
			//only admins can add new users
			else if (P(acl,PERMISSION_OWNERSHIP) || P(acl,PERMISSION_MANAGEUSERS))
			{
				s = con.prepareStatement("INSERT INTO pbrowsedb.membership "
						+ "(groupname, username, access_level) "
						+ "VALUES(?,?, ?);");
				s.setString(1, group);
				s.setString(2, user);
				s.setInt(3, access_level);
				
				s.executeUpdate();
				success = true;
				
				//invalidate
				cm().groupCache.remove(user);
			}
		} 
		catch (SQLException e) 
		{
			//constraint check failure - discard
			success = false;
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
		}
		
		return success;
	}
	
	/**
	 * Removes a user from the mentioned group with the privileges of the 'caller'
	 * @param caller - the user invoking the action
	 * @param user - the user being removed
	 * @param groupname - the group affected
	 * @return success if the caller has user management privilege
	 * @throws SQLException
	 */
	public boolean removeUserFromGroup(String caller, String user, String groupname) throws SQLException
	{
		PreparedStatement s = null;
		boolean success = false;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			int caller_acl = getGroupAccessLevel(caller, groupname);
			
			//only admins can remove users -- or user removes self
			if (P(caller_acl,PERMISSION_OWNERSHIP) || 
					P(caller_acl,PERMISSION_MANAGEUSERS) || 
					caller.equals(user))
			{
				s = con.prepareStatement("DELETE FROM pbrowsedb.membership "
						+ "WHERE username = ? AND groupname = ?;");
				s.setString(1, user);
				s.setString(2, groupname);
				
				s.executeUpdate();
				success = true;
				
				//invalidate
				cm().groupCache.remove(user);
			}
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
			success = false;
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
		}
		
		return success;
	}
	
	/**
	 * Modifies the access level of the 'user' according to the privileges of the 'caller'
	 * @param caller - user invoking the action
	 * @param user - the affected user
	 * @param groupname - the group
	 * @param newlevel - the new access level
	 * @return success if the users privileges were updated
	 * @throws SQLException
	 */
	public boolean updateUserAccessLevel(String caller, String user, String groupname, int newlevel) throws SQLException
	{
		PreparedStatement s = null;
		boolean success = false;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			int caller_acl = getGroupAccessLevel(caller, groupname);
			int user_acl = getGroupAccessLevel(user, groupname);
			
			//cannot set ownership permission
			if ( P(newlevel,PERMISSION_OWNERSHIP) )
			{
				System.out.println("CANNOT SET OWNERSHIP");
				success = false;
			}
			//owner's permissions cannot be changed - by anyone
			else if ( P(user_acl,PERMISSION_OWNERSHIP) )
			{
				System.out.println("CANNOT CHANGE OWNER LEVEL PERMISSIONS");
				success = false;
			}
			//only owners can set this permission
			else if ( P(newlevel,PERMISSION_MANAGEUSERS) && !P(caller_acl,PERMISSION_OWNERSHIP) )
			{
				System.out.println("CANNOT SET MU AS NOT OWNER.");
				success = false;
			}
			//only admins can change user privileges
			else if ( ( P(caller_acl,PERMISSION_OWNERSHIP) || P(caller_acl,PERMISSION_MANAGEUSERS) ) 
					&& !caller.equals(user))
			{
				s = con.prepareStatement("UPDATE pbrowsedb.membership "
						+ "SET access_level = ? WHERE username = ? AND groupname = ?;");
				s.setInt(1, newlevel);
				s.setString(2, user);
				s.setString(3, groupname);
				
				s.executeUpdate();
				success = true;
				
				//invalidate
				cm().groupCache.remove(user);
			}
			else
			{
				System.out.println("NOT OWNER, NOT MU, OR SELF MODIFICATION ATTEMPT.");
				success = false;
			}
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
			success = false;
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
		}
		
		return success;
	}
	
	/**
	 * Adds the file identified by the fileid to the selected group with specified access level
	 * @param caller - the user invoking the request
	 * @param fileid - the file
	 * @param groupname - the affected group
	 * @param access_required - unused parameter
	 * @return success if the file was added
	 * @throws SQLException
	 */
	public boolean addFileACLToGroup(String caller, long fileid, String groupname, int access_required) throws SQLException
	{
		PreparedStatement s0 = null;
		PreparedStatement s1 = null;
		PreparedStatement s2 = null;

		boolean success = false;
		ResultSet r = null;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			int caller_acl = getGroupAccessLevel(caller, groupname);
			
			if (access_required > caller_acl) access_required = caller_acl;
			
			//level 0 users can't do anything to change the group
			if ( P(caller_acl,PERMISSION_OWNERSHIP) || P(caller_acl,PERMISSION_MANAGEFILES) )
			{
				//check ownership of the file
				s0 = con.prepareStatement("SELECT owner, ispublic, indexedby "
						+ "FROM pbrowsedb.data "
						+ "WHERE id = ?");
				s0.setLong(1, fileid);
				
				String owner = null;
				boolean ispublic = false;
				long indexid = 0;
				
				r = s0.executeQuery();
				if (r.next())
				{
					owner = r.getString(1);
					ispublic = r.getBoolean(2);
					indexid = r.getLong(3);
				}
				
				if (owner.equals(caller) || ispublic == true)
				{
					s1 = con.prepareStatement("INSERT INTO pbrowsedb.groupacls "
							+ "(groupname, fileid, require_level) "
							+ "VALUES(?,?,?);");
					s1.setString(1, groupname);
					s1.setLong(2, fileid);
					s1.setLong(3, access_required);
					
					s1.executeUpdate();
					
					//the file has an index, add it too
					if (indexid != 0)
					{
						s2 = con.prepareStatement("INSERT INTO pbrowsedb.groupacls "
								+ "(groupname, fileid, require_level) "
								+ "VALUES(?,?,?);");
						s2.setString(1, groupname);
						s2.setLong(2, indexid);
						s2.setLong(3, access_required);
						s2.executeUpdate();
						
						//invalidate cached index file
						cm().ddCache.remove(indexid);
					}
					
					success = true;
					
					//invalidate cached dd's for affected file(s)
					cm().ddCache.remove(fileid);
				}
			}
		} 
		catch (SQLException e) 
		{
			//might fail due to violation of unique constraint
			//e.g. user tries to add same file twice to same group
			success = false;
		}
		finally 
		{
			if (con != null) con.close();
			if (s0 != null) s0.close();
			if (s1 != null) s1.close();
			if (s2 != null) s2.close();
		}
		
		return success;
	}
	
	/**
	 * Removed the file given by fileid from the access list of the indicated group
	 * @param caller - the user invoking the function
	 * @param fileid - the file
	 * @param groupname - affected group
	 * @return success if the file was removed. The caller must have the ownership, or
	 * file management privilege.
	 * @throws SQLException
	 */
	public boolean deleteFileACLFromGroup(String caller, long fileid, String groupname) throws SQLException
	{
		PreparedStatement s = null;
		boolean success = false;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			int caller_acl = getGroupAccessLevel(caller, groupname);
			
			//uploader or anyone with higher privilege level can delete
			if ( P(caller_acl,PERMISSION_OWNERSHIP) || P(caller_acl,PERMISSION_MANAGEFILES) )
			{
				//remove both the file and it's index from group acls
				s = con.prepareStatement("DELETE FROM pbrowsedb.groupacls "
						+ "WHERE groupname = ? AND "
						+ "( fileid = ? OR "
						+ 	"fileid = ( SELECT indexedby FROM pbrowsedb.data WHERE id = ? ) "
						+ ")");
				s.setString(1, groupname);
				s.setLong(2, fileid);
				s.setLong(3, fileid);
				
				s.executeUpdate();
				success = true;
				
				//invalidate cache
				cm().ddCache.remove(fileid);
			}
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
			success = false;
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
		}
		
		return success;
	}
	
	/**
	 * Get all DDs of files accessible to this group.
	 * @param groupname - the group in question
	 * @return a JsonArray of DataDescriptors for the files accessible by the 
	 * given group
	 * @throws SQLException
	 */
	public JsonArray getGroupFileDDs(String groupname) throws SQLException
	{
		PreparedStatement s = null;
		ResultSet r = null;
		JsonArray result = null;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			s = con.prepareStatement("SELECT d.id, d.owner, g.require_level, d.description, d.path, d.trackname, d.format, d.ispublic, d.indexedby, d.isindex, d.remote "
					+ "FROM pbrowsedb.groupacls g, pbrowsedb.data d "
					+ "WHERE g.groupname = ? AND g.fileid = d.id AND d.isindex = 0");
			s.setString(1, groupname);
	
			result = new JsonArray();
			r = s.executeQuery();
			while (r.next())
			{
				JsonObject jo = new DataDesc(
						r.getLong(1),
						r.getString(2),
						r.getInt(3),
						r.getString(4),
						r.getString(5),
						r.getString(6),
						r.getString(7),
						r.getInt(8), 
						r.getLong(9), 
						r.getBoolean(10), 
						null, r.getBoolean(11))
				.getJson();
				
				result.add(jo);
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
			result = null;
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
			if (r != null) r.close();
		}
		return result;
	}
	
	/**
	 * Returns all the groups the calling user belongs to
	 * @param caller - calling user
	 * @return JsonObject map of group memberships for the given user
	 * @throws SQLException
	 */
	public JsonObject getUserGroups(String caller) throws SQLException
	{
		JsonObject result = null;
		
		HashMap<String,Integer> counts = null;
		HashMap<String,Integer> f_counts = null;
		PreparedStatement s = null;
		ResultSet r = null;
		Connection con = null;
		
		if (cm().groupCache.isKeyInCache(caller))
        {
        	Element e = cm().groupCache.get(caller);
        	if (e != null)
        	{
	        	if (!e.isExpired())
	        	{
	        		result = (JsonObject) e.getObjectValue();
	        	}
        	}
        }

        //the groups are not yet cached, attempt to retrieve them and cache the result
        if (result == null)
        {
        	try {
        		//get group counts
        		counts = this.getGroupUserCounts();
        		f_counts = this.getGroupFileCounts();
        		
        		con = dataSource.getConnection();
        		s = con.prepareStatement("SELECT groupname, access_level "
        				+ "FROM pbrowsedb.membership "
        				+ "WHERE username = ?");
        		s.setString(1, caller);
        		
        		r = s.executeQuery();
        		result = new JsonObject();
        		while (r.next())
        		{
        			JsonObject props = new JsonObject();
        			props.addProperty("access_level", r.getInt(2));
        			props.addProperty("count", counts.get(r.getString(1)));
        			props.addProperty("file_count", f_counts.get(r.getString(1)));
        			
        			result.add(r.getString(1), props);
        		}
        		cm().groupCache.put(new Element(caller,result));
        	}
        	catch (Exception e)
        	{
        		e.printStackTrace();
        		result = null;
        	}
        	finally 
        	{
        		if (con != null) con.close();
        		if (s != null) s.close();
        		if (r != null) r.close();
        	}
        }
		
		return result;
	}
	
	/**
	 * Gets a mapping of all the groups with access to the given file
	 * @param fileid - the file in question
	 * @return JsonObject mapping of groups
	 * @throws SQLException
	 */
	public JsonObject getFileGroups(long fileid) throws SQLException
	{
		JsonObject result = new JsonObject();
		PreparedStatement s = null;
		ResultSet r = null;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			s = con.prepareStatement("SELECT groupname, require_level "
					+ "FROM pbrowsedb.groupacls "
					+ "WHERE fileid = ?");
			s.setLong(1, fileid);
	
			r = s.executeQuery();
			while (r.next())
			{
				result.addProperty(r.getString(1),r.getInt(2));
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
			result = null;
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
			if (r != null) r.close();
		}
		return result;
	}
	
	/**
	 * Gets the number of files belonging to every group
	 * @return Map of group to the number of files in that group
	 * @throws SQLException
	 */
	public HashMap<String,Integer> getGroupFileCounts() throws SQLException
	{
		HashMap<String,Integer> result = new HashMap<String,Integer>();
		PreparedStatement s = null;
		ResultSet r = null;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			s = con.prepareStatement("SELECT g.groupname, count(groupname) "
					+ "FROM pbrowsedb.groupacls g, pbrowsedb.data d "
					+ "WHERE g.fileid = d.id AND d.isindex = 0 "
					+ "GROUP BY g.groupname;");
			r = s.executeQuery();
			
			while (r.next())
			{
				result.put(r.getString(1),r.getInt(2));
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
			result = null;
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
			if (r != null) r.close();
		}
		return result;
	}
	
	/**
	 * Get the number of users belonging to every group
	 * @return a map of group to number of users in that group
	 * @throws SQLException
	 */
	public HashMap<String,Integer> getGroupUserCounts() throws SQLException
	{
		HashMap<String,Integer> result = new HashMap<String,Integer>();
		PreparedStatement s = null;
		ResultSet r = null;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			s = con.prepareStatement("SELECT groupname, count(groupname) "
					+ "FROM pbrowsedb.membership GROUP BY groupname;");
			r = s.executeQuery();
			
			while (r.next())
			{
				result.put(r.getString(1),r.getInt(2));
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
			result = null;
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
			if (r != null) r.close();
		}
		return result;
	}
	
	/**
	 * Gets a list of all the users belonging the given group along with their privilege levels
	 * @param groupname - the group in question
	 * @return JsonArray of users with their access_level
	 * @throws SQLException
	 */
	public JsonArray getGroupUsers(String groupname) throws SQLException
	{
		JsonArray result = new JsonArray();
		
		PreparedStatement s = null;
		ResultSet r = null;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			s = con.prepareStatement("SELECT username, access_level "
					+ "FROM pbrowsedb.membership "
					+ "WHERE groupname = ?");
			s.setString(1, groupname);
			r = s.executeQuery();
			
			while (r.next())
			{
				System.out.println(r.getString(1));
				
				JsonObject prop = new JsonObject();
				prop.addProperty("username", r.getString(1));
				prop.addProperty("access_level", r.getInt(2));
				
				result.add(prop);
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
			result = null;
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
			if (r != null) r.close();
		}

		System.out.println(result.toString());
		return result;
	}
	
	/**
	 * Stores the uauth key of a logged in user in the database such that their session can
	 * be resumed without needing to login again
	 * @param username - the affected user
	 * @param uauth - the uauth token generated at login
	 * @return success if the token was stored
	 * @throws SQLException
	 */
	public boolean storeSessionKey(String username, String uauth) throws SQLException
	{
		if (username == null || uauth == null)
		{
			return false;
		}
		
		PreparedStatement s = null;
		Connection con = null;
		boolean result = false;
		
		try {
			con = dataSource.getConnection();
			s = con.prepareStatement("UPDATE pbrowsedb.users SET sessionkey = ? "
				+ "WHERE username = ? ");
			s.setString(1, uauth);
			s.setString(2, username);
			s.executeUpdate();
			result = true;
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
		}
		
		return result;
	}
	
	/**
	 * Verifies the integrity of a user supplied uauth when attempting to resume a previously
	 * interrupted session
	 * @param caller - the invoking user
	 * @param uauth - the authentication token
	 * @return login information as if the user were loggin in normally
	 * @throws SQLException
	 */
	public JsonObject checkSessionKey(String caller, String uauth) throws SQLException
	{
		if (caller == null)
		{
			return null;
		}
		
		PreparedStatement s = null;
		Connection con = null;
		ResultSet r = null;
		JsonObject result = null;
		
		try {
			con = dataSource.getConnection();
			s = con.prepareStatement("SELECT username, nick "
					+ "FROM pbrowsedb.users "
					+ "WHERE username = ? AND sessionkey = ?");
			s.setString(1, caller);
			s.setString(2, uauth);
			r = s.executeQuery();
			
			if (r.next())
			{
				result = new JsonObject();
				result.addProperty("username", r.getString(1));
				result.addProperty("nick", r.getString(2));
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
			result = null;
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
			if (r != null) r.close();
		}
		return result;
	}
	
	/**
	 * Gets a map of the emails of all the group members of the given group
	 * @param caller - the invoking user
	 * @param group - the group in question
	 * @return A map of username to email if successful. Requires the invoking user to have
	 * the group ownership privilege.
	 * @throws SQLException
	 */
	public HashMap<String,String> getGroupEmails(String caller, String group) throws SQLException
	{
		HashMap<String,String> emails = null;
		
		PreparedStatement s = null;
		Connection con = null;
		ResultSet r = null;
		
		try 
		{
			int caller_acl = getGroupAccessLevel(caller, group);
		
			//only owner can get emails
			if (P(caller_acl,PERMISSION_OWNERSHIP))
			{
				con = dataSource.getConnection();
				s = con.prepareStatement("SELECT u.username, u.email "
						+ "FROM pbrowsedb.membership m, pbrowsedb.users u "
						+ "WHERE u.username = m.username AND m.groupname = ?;");
				s.setString(1, group);
				r = s.executeQuery();
				
				emails = new HashMap<String,String>();
				while (r.next())
				{
					emails.put(r.getString(1), r.getString(2));
				}
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
			if (r != null) r.close();
		}
		
		return emails;
	}
	
	/**
	 * Changes the password of the given user to the newly provided one
	 * @param caller - the affected user
	 * @param old - the old password, for verification
	 * @param newpass - the password to change to
	 * @return success if the password was changed
	 * @throws SQLException
	 */
	public boolean changeUserPassword(String caller, String old, String newpass) throws SQLException
	{
		if (caller == null || old == null || newpass == null)
		{
			return false;
		}
		
		PreparedStatement s = null;
		PreparedStatement s2 = null;
		Connection con = null;
		boolean result = false;
		ResultSet r = null;
		
		try {
			con = dataSource.getConnection();

			s = con.prepareStatement("SELECT password "
					+ "FROM pbrowsedb.users "
					+ "WHERE username = ?;");
			s.setString(1, caller);
			
			String _old = null;
			r = s.executeQuery();
			if (r.next())
			{
				_old = r.getString(1);
			}
			
			if (_old.equals(old))
			{
				s2 = con.prepareStatement("UPDATE pbrowsedb.users SET password = ? "
						+ "WHERE username = ? ");
				s2.setString(1, newpass);
				s2.setString(2, caller);
				s2.executeUpdate();
				result = true;
			}
			else result = false;
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
			if (s2 != null) s2.close();
			if (r != null) r.close();
		}
		
		return result;
	}
	
	/**
	 * Forcibly changes the user's password to the new one
	 * @param caller - affected user
	 * @param newpass - new hashed password
	 * @return success
	 * @throws SQLException
	 */
	public boolean resetUserPassword(String caller, String newpass) throws SQLException
	{
		if (caller == null || newpass == null)
		{
			return false;
		}
		
		PreparedStatement s = null;
		Connection con = null;
		boolean result = false;
		ResultSet r = null;
		
		try {
			con = dataSource.getConnection();

			s = con.prepareStatement("UPDATE pbrowsedb.users SET password = ? "
					+ "WHERE username = ? ");
			s.setString(1, newpass);
			s.setString(2, caller);
			s.executeUpdate();
			result = true;
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
			if (r != null) r.close();
		}
		
		return result;
	}

	/**
	 * Changes the email address of the requesting user
	 * @param username - invoking user
	 * @param email - new email address
	 * @return success if the request was processed. The user is required to access the
	 * token sent to the new email address in order to complete the action.
	 * @throws SQLException
	 */
	public boolean changeUserEmail(String username, String email) throws SQLException 
	{
		if (username == null || email == null)
		{
			return false;
		}
		
		PreparedStatement s = null;
		Connection con = null;
		boolean result = false;
		
		try {
			con = dataSource.getConnection();

			s = con.prepareStatement("UPDATE pbrowsedb.users SET email = ? "
					+ "WHERE username = ? ");
			s.setString(1, email);
			s.setString(2, username);
			s.executeUpdate();
			result = true;
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
		}
		
		return result;
	}
	
	/**
	 * Sets an authentication token request for the given user, performing the specified task
	 * @param username - the user
	 * @param token - the authentication token
	 * @param function - the operation which this token performs
	 * @param params - the parameters required to perform the operation
	 * @return success if the token was set
	 * @throws SQLException
	 */
	public boolean setToken(String username, String token, int function, String params) throws SQLException
	{
		if (username == null || token == null)
		{
			return false;
		}
		
		if (token.length() != 32)
		{
			return false;
		}
		
		PreparedStatement s = null;
		Connection con = null;
		boolean result = false;
		
		try {
			con = dataSource.getConnection();

			s = con.prepareStatement("INSERT INTO pbrowsedb.tokens (username,token,funct,params) "
					+ "VALUES(?,?,?,?);");
			s.setString(1, username);
			s.setString(2, token);
			s.setInt(3, function);
			s.setString(4, params);
			s.executeUpdate();
			result = true;
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
		}
		
		return result;
	}
	
	/**
	 * Retrieves a token with the given authentication hash
	 * @param hash - the token
	 * @return Map of token properties required to perform tokenised action
	 * @throws SQLException
	 */
	public HashMap<String,Object> getToken(String hash) throws SQLException
	{
		HashMap<String,Object> token = null;
		
		if (hash == null)
		{
			return null;
		}
		
		PreparedStatement s = null;
		PreparedStatement s2 = null;
		Connection con = null;
		ResultSet r = null;
		
		try {
			con = dataSource.getConnection();

			s = con.prepareStatement("SELECT id, username, funct, params "
					+ "FROM pbrowsedb.tokens "
					+ "WHERE token = ?;");
			s.setString(1, hash);
			
			long id = -1;
			
			r = s.executeQuery();
			if (r.next())
			{
				token = new HashMap<String,Object>();
				token.put("username", r.getString(2));
				token.put("funct", r.getInt(3));
				token.put("params", r.getString(4));
				id = r.getLong(1);
			}
			
			if (id != -1)
			{
				//delete the token
				s2 = con.prepareStatement("DELETE FROM pbrowsedb.tokens "
						+ "WHERE token = ?;");
				s2.setString(1, hash);
				s2.executeUpdate();
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
			token = null;
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
			if (s2 != null) s2.close();
			if (r != null) r.close();
		}
		return token;
	}

	/**
	 * Changes the nickname of the calling user to the newly specified one
	 * @param username - calling user
	 * @param nickname - new nickname
	 * @return success if the nickname was changed
	 * @throws SQLException
	 */
	public boolean changeUserNickname(String username, String nickname) throws SQLException
	{
		if (username == null || nickname == null)
		{
			return false;
		}
		
		PreparedStatement s = null;
		Connection con = null;
		boolean result = false;
		
		try {
			con = dataSource.getConnection();

			s = con.prepareStatement("UPDATE pbrowsedb.users SET nick = ? "
					+ "WHERE username = ? ");
			s.setString(1, nickname);
			s.setString(2, username);
			s.executeUpdate();
			result = true;
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
		}
		
		return result;
	}
	
	/**
	 * Sets a user account as being activated - required for login. This function is invoked by
	 * the registration token task.
	 * @param username - the affected user
	 * @return success if the account was activated
	 * @throws SQLException
	 */
	public boolean activateUserAccount(String username) throws SQLException
	{
		if (username == null)
		{
			return false;
		}
		
		PreparedStatement s = null;
		Connection con = null;
		boolean result = false;
		
		try {
			con = dataSource.getConnection();

			s = con.prepareStatement("UPDATE pbrowsedb.users SET activated = true "
					+ "WHERE username = ? ");
			s.setString(1, username);
			s.executeUpdate();
			result = true;
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
		}
		
		return result;
	}
	
	/**
	 * Checks the existence of a particular user from either its email address or username
	 * If the user exists, returns both the email and username of the user
	 * @param usermail
	 * @return
	 * @throws SQLException
	 */
	public HashMap<String,String> lookupUser(String usermail) throws SQLException
	{
		HashMap<String,String> resultset = null;
		
		PreparedStatement s = null;
		Connection con = null;
		ResultSet r = null;
		
		try {
			con = dataSource.getConnection();

			s = con.prepareStatement("SELECT username, email "
					+ "FROM pbrowsedb.users "
					+ "WHERE username = ? OR email = ?;");
			s.setString(1, usermail);
			s.setString(2, usermail);
			
			r = s.executeQuery();
			if (r.next())
			{
				resultset = new HashMap<String,String>();
				resultset.put("username", r.getString(1));
				resultset.put("email", r.getString(2));
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
			if (r != null) r.close();
		}
		
		return resultset;
	}
	
	/**
	 * Stores the passed in user config string, as a stringified JSON parameter.
	 * @param username - the calling user
	 * @param config - the config string
	 * @return
	 * @throws SQLException
	 */
	public boolean updateStoredConfig(String username, String config) throws SQLException
	{
		if (username == null)
		{
			System.err.println("No user provided.");
			return false;
		}
		
		PreparedStatement s = null;
		boolean success = true;
		Connection con = null;
		
		try {
			con = dataSource.getConnection();
			
			s = con.prepareStatement("UPDATE pbrowsedb.users SET saved_config = ? "
					+ "WHERE username = ? ");
			s.setString(1, config);
			s.setString(2, username);
			s.executeUpdate();
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
			success = false;
		}
		finally 
		{
			if (con != null) con.close();
			if (s != null) s.close();
		}
		return success;
	}
}
