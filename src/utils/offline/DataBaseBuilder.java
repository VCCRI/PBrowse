package utils.offline;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;

import com.database.mysql.DBA;
import com.database.mysql.DataDesc;

import static com.database.mysql.DBA.dba;

/**
 * Parse all files in the provided directory directory and completely rebuild all data entries - according
 * to the provided metadata file.
 * @author root
 *
 */
public class DataBaseBuilder 
{
	public static final String path = "/opt/pbrowse/sessions/root/";
	public static final String metapath = "/root/aspdata/Table4PBrowse.csv";
	private ArrayList<File> allfiles;
	
	public DataBaseBuilder()
	{
		allfiles = new ArrayList<File>();
	}
	
	private void explore(File f)
	{
		if (!f.isDirectory())
		{
			allfiles.add(f);
			return;
		}
		
		for (File child : f.listFiles())
		{
			this.explore(child);
		}
	}
	
	public void run(OutputStreamWriter out) throws SQLException, IOException 
	{
		//quick recursive search for all files in the sub-directories
		File root = new File(path);
		for (File child : root.listFiles())
		{
			explore(child);
		}
		
		//read metadata file and build data structure
		BufferedReader br = new BufferedReader(new FileReader(metapath));
		String line = null; //ignore first
		
		HashMap<String,DataDesc> map = new HashMap<String,DataDesc>();
		
		while ((line = br.readLine()) != null)
		{
			if (line.length() < 1)
			{
				continue;
			}
			
			String data[] = line.split("\t");
			
			String trackname = data[0];
			String desc = data[1];
			String path = data[2];
			
			String p[] = path.split("/");
			String filename = p[p.length-1];
			out.write(filename+"\n");
			
			map.put(filename, new DataDesc(0,"",0,desc,"",trackname,"BigWig",DBA.PRIVACY_PUBLIC,0,false,null,false));
		}
		br.close();
		
		ArrayList<String> residual = new ArrayList<String>();
		
		//build dd for each file
		for (File f : allfiles)
		{
			String p = f.getAbsolutePath();
			String parts[] = p.split("/");
			
			String username = parts[4];
			out.write("username: "+username+"\n");
			
			String fpath = "";
			for (int i=5; i<parts.length; i++)
			{
				fpath += "/"+parts[i];
			}
			String filename = parts[parts.length-1];
			
			DataDesc dd = map.get(filename);
			if (dd == null)
			{
				residual.add(filename);
				continue;
			}
			
			out.write("filepath: "+fpath+"\n");
			out.write("filename: "+filename+"\n");
			
			String fmt[] = parts[parts.length-1].split("\\.");
			String format = fmt[fmt.length-1];
			
			switch (format.toLowerCase())
			{
			case "bw":
				format = "BigWig";
				break;
			case "bb":
				format = "BigBed";
				break;
			case "bigbed":
				format = "BigBed";
				break;
			case "bigwig":
				format = "BigWig";
				break;
			}
			
			out.write("format: "+format+"\n");
			
			if (dba().lookupUser(username) == null)
			{
				//create a new user for the purpose
				dba().doRegister(username, username, "root@pbrowse.com", "root");
			}
			
			dba().registerDataFile(username, dd.description, dd.trackname, fpath, format, DBA.PRIVACY_PUBLIC, null, false, false);
		}
		
		out.write("Remaining: "+"\n");
		for (String s : residual)
		{
			out.write(s+"\n");
		}
	}

}
