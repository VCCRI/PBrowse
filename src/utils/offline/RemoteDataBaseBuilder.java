package utils.offline;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.sql.SQLException;
import java.util.HashMap;

import com.database.mysql.DBA;
import com.google.gson.JsonObject;

import static com.database.mysql.DBA.dba;

/**
 * Parse all files in the provided directory directory and completely rebuild all data entries - according
 * to the provided metadata file.
 * @author root
 *
 */
public class RemoteDataBaseBuilder 
{
	public static final String metapath = "/opt/pbrowse/trackDb.txt";
	
	public RemoteDataBaseBuilder()
	{
		
	}
	
	public void run(OutputStreamWriter out) throws SQLException, IOException 
	{
		//read metadata file and build data structure
		String line = null; //ignore first
		
		//get first level no tab tracks
		//track -> trackprops-map
		HashMap<String,HashMap<String,String>> level0 = new HashMap<String,HashMap<String,String>>();
		
		BufferedReader br = new BufferedReader(new FileReader(metapath));
		String c_track = null;
		while ((line = br.readLine()) != null)
		{
			if (line.length() < 1 || line.startsWith("\t") )
			{
				c_track = null;
				continue;
			}
			
			String data[] = line.trim().split(" ");
			if (data[0].equals("track"))
			{
				c_track = data[1];
			}
			
			if (c_track != null)
			{
				if (!level0.containsKey(c_track))
				{
					level0.put(c_track, new HashMap<String,String>());
				}
				
				if (data.length == 1)
				{
					continue;
				}
				
				HashMap<String,String> props = level0.get(c_track);
				
				String combined = "";
				for (int i=1; i<data.length; i++)
				{
					combined += data[i]+" ";
				}
				props.put(data[0], combined.trim());
			}
		}
		br.close();
		
		//level1 tracks
		HashMap<String,HashMap<String,String>> level1 = new HashMap<String,HashMap<String,String>>();
		
		br = new BufferedReader(new FileReader(metapath));
		c_track = null;
		while ((line = br.readLine()) != null)
		{
			if (!line.startsWith("\t") || line.startsWith("\t\t"))
			{
				c_track = null;
				continue;
			}
			
			String data[] = line.trim().split(" ");
			if (data[0].equals("track"))
			{
				c_track = data[1];
			}
			
			if (c_track != null)
			{
				if (!level1.containsKey(c_track))
				{
					level1.put(c_track, new HashMap<String,String>());
				}
				
				if (data.length == 1)
				{
					continue;
				}
				
				HashMap<String,String> props = level1.get(c_track);
				
				String combined = "";
				for (int i=1; i<data.length; i++)
				{
					combined += data[i]+" ";
				}
				props.put(data[0], combined.trim());
			}
		}
		br.close();
		
		//level1 tracks
		HashMap<String,HashMap<String,String>> level2 = new HashMap<String,HashMap<String,String>>();
		
		br = new BufferedReader(new FileReader(metapath));
		c_track = null;
		while ((line = br.readLine()) != null)
		{
			if (!line.startsWith("\t\t"))
			{
				c_track = null;
				continue;
			}
			
			String data[] = line.trim().split(" ");
			if (data[0].equals("track"))
			{
				c_track = data[1];
			}
			
			if (c_track != null)
			{
				if (!level2.containsKey(c_track))
				{
					level2.put(c_track, new HashMap<String,String>());
				}
				
				if (data.length == 1)
				{
					continue;
				}
				
				HashMap<String,String> props = level2.get(c_track);
				
				String combined = "";
				for (int i=1; i<data.length; i++)
				{
					combined += data[i]+" ";
				}
				props.put(data[0], combined.trim());
			}
		}
		br.close();
		
//		for (String s : level0.keySet())
//		{
//			out.write(s+"\n");
//		}
		
		for (String s : level1.keySet())
		{
			out.write(s+" "+level1.get(s).get("parent")+"\n");
		}
		
		if (dba().lookupUser("ENCODE") == null)
		{
			//create a new user for the purpose
			dba().doRegister("ENCODE", "ENCODE", "root@pbrowse.com", "ENCODE");
		}
		
		for (String s : level2.keySet())
		{
			try
			{
				HashMap<String,String> props = level2.get(s);
				if (props.get("parent") == null)
				{
					props.put("parent", "");
				}
				String myparent = props.get("parent").split(" ")[0];
				out.write(myparent+"\n");
				
				String parentparent = null;
				if (myparent.length() == 0 )
				{
					parentparent = "";
				}
				else parentparent = level1.get(myparent).get("parent");
				
				String desc = props.get("longLabel")+" "+props.get("subGroups")+" "+props.get("metadata");
				String trackname = props.get("shortLabel");
				String path = props.get("bigDataUrl");
				
				//build format meta
				JsonObject format = new JsonObject();
				format.addProperty("studyid", parentparent+"/"+myparent);
				
				if (props.get("type") == null)
				{
					format.addProperty("data_format", "");
				}
				else format.addProperty("data_format", props.get("type").split(" ")[0].toLowerCase());
				
				format.addProperty("filename", props.get("track"));
				format.addProperty("pennant", "http://genome.ucsc.edu/images/encodeThumbnail.jpg");
				
				out.write(desc+"\n"+trackname+"\n"+path+"\n"+format.toString()+"\n"+"\n");
				dba().registerDataFile("ENCODE", desc, trackname, path, format.toString(), DBA.PRIVACY_PUBLIC, null, false, true);
				
			} 
			catch (Exception e)
			{
				e.printStackTrace();break;
			}
			
		}
		out.write(level2.keySet().size()+"\n");
		
		
	}

}
