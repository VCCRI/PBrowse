package com.backend.filesystem;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.sql.SQLException;

import javax.servlet.ServletException;
import javax.servlet.annotation.MultipartConfig;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.Cookie;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.servlet.http.Part;

import com.backend.collab.IUser;
import com.backend.collab.PSession;
import com.backend.collab.SessionManager;
import com.database.mysql.DBA;
import com.google.gson.JsonObject;

import static com.database.mysql.DBA.dba;

/**
 * Servlet implementation class UploadEndPoint
 */
@WebServlet("/uploadData")
@MultipartConfig
public class UploadEndpoint extends HttpServlet 
{
	private static final long serialVersionUID = 1L;
	private final String basedir = "/opt/pbrowse/";
	
    public UploadEndpoint() {}
	protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {}

	private String getFileName(Part fp)
	{
		return fp.getHeader("content-disposition").split("=")[2].replaceAll("\"", "");
	}
	
	private final int NOINDEX = 8;
	private final int ZIPPED = 4;
	private final int INDEXFILE = 2;
	private final int DATAFILE = 1;

	private final String [] filetypes = { "bigbed", "bb", "bigwig", "bw", "2bit", "bam", "vcf", "bed", "wig" };
	private final String [] indexes = { "bai", "tbi" };
	
	/**
	 * Returns the last extension of a file, indicating its most probable type,
	 * ignores gz ending
	 * @param s
	 * @return
	 */
	private String getFileExt(String s)
	{
		String exts[] = s.split("\\.");
		for (int i=exts.length-1; i>=0; i--)
	    {
	    	for (int j=0; j<filetypes.length; j++)
	    	{
	    		if (exts[i].equalsIgnoreCase(filetypes[j]))
	    		{
	    			return exts[i];
	    		}
	    	}
	    	for (int j=0; j<indexes.length; j++)
	    	{
	    		if (exts[i].equalsIgnoreCase(indexes[j]))
	    		{
	    			return exts[i];
	    		}
	    	}
	    }
		return "";
	}
	
	/**
	 * Returns a set of flags indicating the type of file being processed - i.e. Index or Data
	 * @param s
	 * @return
	 */
	private int getFileType(String s)
	{
		if (s.length() == 0 || s == null)
		{
			return NOINDEX;
		}
		String exts[] = s.split("\\.");
		
	    int zipped = 0, validFile = 0, isindex = 0;
	    
	    for (int i=exts.length-1; i>=0; i--)
	    {
	    	for (int j=0; j<filetypes.length; j++)
	    	{
	    		if (exts[i].equalsIgnoreCase(filetypes[j]))
	    		{
	    			validFile = DATAFILE;
	    			break;
	    		}
	    	}
	    	for (int j=0; j<indexes.length; j++)
	    	{
	    		if (exts[i].equalsIgnoreCase(indexes[j]))
	    		{
	    			isindex = INDEXFILE;
	    			break;
	    		}
	    	}
	    	
	    	if (exts[i].equalsIgnoreCase("gz"))
    		{
    			zipped = ZIPPED;
    		}
	    }
	    
	    return (zipped | validFile | isindex);
	}

	//Don't unpack the file if it is gzipped
	private boolean writeStandard(File f, Part p)
	{
	    try (InputStream input = p.getInputStream()) {
	    	Files.copy(input, f.toPath(), StandardCopyOption.REPLACE_EXISTING);
	    }
	    catch (Exception e)
	    {
	    	return false;
	    }
	    return true;
	}
	
	/**
	 * Handles uploading of genomic data for sharing with other users
	 */
	protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException 
	{
		JsonObject resp = new JsonObject();
		OutputStreamWriter os = new OutputStreamWriter(response.getOutputStream());
		
	    Cookie[] cookies = request.getCookies();
	    if (cookies == null)
	    {
	    	resp.addProperty("status", "failed");
	    	resp.addProperty("error", "No cookie data.");
	    	
	    	os.write(resp.toString());
	    	os.close();
	    	return;
	    }
	    
	    //get authentication data from cookies
	    String uauth = null, username = null;
	    for (Cookie c : cookies)
	    {
	    	if (c.getName().equals("uauth"))
	    	{
	    		uauth = c.getValue();
	    	}
	    	else if (c.getName().equals("username"))
	    	{
	    		username = c.getValue();
	    	}
	    }
	    //if a field is missing, assume bad authentication
	    if (uauth == null || username == null)
	    {
	    	resp.addProperty("status", "failed");
	    	resp.addProperty("error", "Bad cookie data.");
	    	
	    	os.write(resp.toString());
	    	os.close();
	    	return;
	    }
	    
	    SessionManager sm = SessionManager.getISM();
	    IUser u = sm.getUserFromName(username);
	    
	    //check provided data matches the entry stored by the server
	    if (!u.getUAuth().equals(uauth) || !u.username.equals(username))
	    {
	    	resp.addProperty("status", "failed");
	    	resp.addProperty("error", "Bad authentication.");
	    	
	    	os.write(resp.toString());
	    	os.close();
	    	return;
	    }

	    //build a path for the new file based on:
	    //{basedir}/sessions/{username}/{studyname}/{filename}
	    String study = request.getParameter("studyid");
	    String desc = request.getParameter("description");
	    String trackname = request.getParameter("trackname");
	    String tag = request.getParameter("genometag");
	    
	    //prefix tag to description
	    desc = tag+(desc.length()>0?": "+desc:"");
	    
	    int ispublic = Integer.parseInt(request.getParameter("ispublic"));
	    
	    //ensure all parameters are provided
	    if (study == null || desc == null || trackname == null)
	    {
	    	return;
	    }

	    //get the actual file
	    Part data_file = request.getPart("data_file"); // Retrieves <input type="file" name="file">
	    Part index_file = request.getPart("index_file"); // Retrieves <input type="file" name="file">

	    //get filename from header info
	    String fileName = getFileName(data_file);
	    String indexFileName = getFileName(index_file);
	    
	    System.out.println("FileName: "+fileName+" index: "+indexFileName);
	    
	    //make directory structure
	    File uploads = new File(this.basedir);
	    File dir_root = new File(uploads, "sessions/"+u.username+"/"+study+"/");
	    try {
	    	dir_root.mkdirs();
	    } catch (Exception e) {
	    	resp.addProperty("status", "failed");
	    	resp.addProperty("error", "Failed to create directories.");
	    	
	    	os.write(resp.toString());
	    	os.close();
	    	return;
	    }
	    
	    //create files for index and data
	    File file = new File(dir_root, fileName);
	    File indexfile = new File(dir_root, indexFileName);

	    //check the data file extensions are acceptable
	    int dr = getFileType(fileName);
	    String dext = getFileExt(fileName);
	    if ((dr & DATAFILE) == 0)
	    {
	    	resp.addProperty("status", "failed");
	    	resp.addProperty("error", "Not a recognised genomic data file! Aborting upload.");
	    	os.write(resp.toString());
	    	os.close();
	    	return;
	    }
	    
	    //Check the index has acceptable format if it is required
	    int ir = getFileType(indexFileName);
	    String iext = getFileExt(indexFileName);
	    if ((ir & INDEXFILE) == 0 && (ir & NOINDEX) == 0)
	    {
	    	resp.addProperty("status", "failed");
	    	resp.addProperty("error", "Not a recognised index file! Aborting upload.");
	    	os.write(resp.toString());
	    	os.close();
	    	return;
	    }
	    
	    //check for a valid combination of extensions
	    boolean goodindex = true, index_required = false;
	    if (dext.equals("bam"))
	    {
	    	if (!iext.equals("bai"))
	    	{
	    		goodindex = false;
	    	}
	    	index_required = true;
	    }
	    if (dext.equals("bed") || dext.equals("vcf"))
	    {
	    	if (!iext.equals("tbi"))
	    	{
	    		goodindex = false;
	    	}
	    	index_required = true;
	    }
	    if (!goodindex)
	    {
	    	resp.addProperty("status", "failed");
	    	resp.addProperty("error", "Wrong index file for data file.");
	    	
	    	os.write(resp.toString());
	    	os.close();
	    	return;
	    }
	    
	    //replace bb and bw extensions with the standard BigBed and BigWig for uniformity
	    if (fileName.contains(".bb"))
	    	fileName.replaceAll("\\.bb", "\\.BigBed");
	    if (fileName.contains(".bw"))
	    	fileName.replaceAll("\\.bw", "\\.BigWig");
	    
	    //register the path of the new file in the database
	    String path = "/"+study+"/"+fileName;
	    boolean registered = true;
	    Integer index_id = null;
	    Integer fileid = null;
	    try {
		    if (index_required)
		    {
		    	String indexpath = "/"+study+"/"+indexFileName;
				if ((index_id = dba().registerDataFile(username, "", "", indexpath, iext, ispublic, null, true, false)) == -1)
				{
					registered = false;
				}
		    }
		    //register the data file indexed by the already registered index file - if required
			if ((fileid = dba().registerDataFile(username, desc, trackname, path, dext, ispublic, index_id, false, false)) == -1)
			{
				registered = false;
			}
		} catch (SQLException e) {
			e.printStackTrace();
			registered = false;
		}
	    
	    if (!registered)
	    {
	    	resp.addProperty("status", "failed");
	    	resp.addProperty("error", "Database could not register file(s).");
	    	
	    	os.write(resp.toString());
	    	os.close();
	    	return;
	    }
	    
	    //actually write the files to disk
	    boolean standard_write = true;
    	if (!writeStandard(file, data_file))
    		standard_write = false;
    	
	    if (index_required)
	    {
	    	if (!writeStandard(indexfile, index_file))
	    		standard_write = false;
	    }
	    if (!standard_write)
	    {
	    	System.out.println("Failed to write file(s)!");
	    	resp.addProperty("status", "failed");
	    	resp.addProperty("error", "Could not write the file(s).");
	    }
	    else 
	    { 
	    	resp.addProperty("status", "success");

	    	//for notifying other users
	    	resp.addProperty("id", fileid);
	    	resp.addProperty("owner", u.username);
	    	resp.addProperty("description", desc);
	    	resp.addProperty("path", path);
	    	resp.addProperty("format", dext);
	    	resp.addProperty("trackname", trackname);
	    	resp.addProperty("indexedby", index_id);
	    	resp.addProperty("ispublic", ispublic);
	    	resp.addProperty("isindex", false);
	    	resp.addProperty("remote", false);
	    	if (ispublic == DBA.PRIVACY_PUBLIC)
	    	{
	    		resp.addProperty("rtype", "new-public-file");
	    		sm.broadcastMessage(resp.toString());
	    	}
    		//inform session users that a new file was uploaded
    		resp.addProperty("rtype", "new-session-file");
    		PSession ps = sm.getUserSession(u);
    		if (ps != null && ispublic == DBA.PRIVACY_PRIVATE)
    		{
    			ps.updateAllUsers(resp.toString());
    		}
    		else u.sendMessage(resp.toString());
    	}

	    System.out.println(resp.toString());
	    
	    os.write(resp.toString());
	    os.close();
	}
}
