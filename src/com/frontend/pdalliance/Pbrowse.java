package com.frontend.pdalliance;

import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;

import javax.servlet.ServletConfig;
import javax.servlet.ServletException;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import com.database.mysql.DBA;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import static com.database.mysql.DBA.dba;
import static com.backend.collab.SessionManager.getISM;

/**
 * Serves the Pbrowse-dalliance html to the end user
 */
@WebServlet("")
public class Pbrowse extends HttpServlet {
	private static final long serialVersionUID = 1L;
       
    public Pbrowse() 
    {
        super();
    }
    
    /**
     * Initialise the servlet, along with all several important singleton manager instances
     */
    public void init(ServletConfig config) 
    {
    	//init cache manager
    	com.backend.collab.CacheManager.cm();
    	
    	//init session manager
    	com.backend.collab.SessionManager.getISM();
    }

    /**
     * Get request handler for pbrowse
     */
    protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException 
	{
    	String token = request.getParameter("token");
    	String testmode = request.getParameter("testing");

    	if (testmode != null)
    	{
    		request.setAttribute("testmode", true);
    		request.getRequestDispatcher("index.jsp").forward(request, response);
    		return;
    	}
    	
		/**
		 * Process requests with the token parameter differently
		 */
		if (token != null)
		{
			HashMap<String,Object> f = null;
			try {
				f = dba().getToken(token);
				if (f != null)
				{
					String username = (String) f.get("username");
					Integer funct 	= (Integer) f.get("funct");
					
					String params 	= (String) f.get("params");
					JsonObject p = new JsonParser().parse(params).getAsJsonObject();
					
					boolean success = false;
					
					switch (funct)
					{
					/**
					 * Changes the email of the specified user
					 */
					case DBA.FUNCTION_CHANGEEMAIL:
						request.setAttribute("type", "change-email");
						success = dba().changeUserEmail(username, p.get("email").getAsString());
						if (success)
						{
							request.setAttribute("status", "success");
							request.setAttribute("error", "The email address has been changed successfully.");
						}
						break;
					/**
					 * Activates the specified user's account
					 */
					case DBA.FUNCTION_REGISTERACCOUNT:
						success = dba().activateUserAccount(username);
						if (success)
						{
							request.setAttribute("status", "success");
							request.setAttribute("error", "Your account has been activated.");
						}
						break;
						
					/**
					 * Activates the specified user's account
					 */
					case DBA.FUNCTION_RESETPASSWORD:
						String newpass = p.get("newpassword").getAsString();
						success = dba().resetUserPassword(username, getISM().sha512(newpass));
						if (success)
						{
							request.setAttribute("status", "success");
							request.setAttribute("error", "Your password has been changed to: '"+newpass+"'. "
									+ "To change it, visit your profile page when logged in.");
						}
						break;
						
					/**
					 * Unrecognised or not yet implemented
					 */
					default:
						request.setAttribute("status", "fail");
						request.setAttribute("error", "The token is invalid, or has been used already.");
						break;
					}
					
				}
				else
				{
					request.setAttribute("status", "fail");
					request.setAttribute("error", "The token is invalid, or has been used already.");
				}
			} catch (SQLException e) {
				e.printStackTrace();
			}
			
			request.getRequestDispatcher("token.jsp").forward(request, response);
			return;
		}
		
		/**
		 * Direct to standard pbrowse interface
		 */
		request.getRequestDispatcher("index.jsp").forward(request, response);
	}

    /**
     * same as get
     */
	protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException 
	{
		doGet(request, response);
	}

}
