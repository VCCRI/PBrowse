package com.frontend.pdalliance;

import java.io.IOException;
import java.io.OutputStreamWriter;
import java.sql.SQLException;

import javax.servlet.ServletException;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import utils.offline.DataBaseBuilder;
import utils.offline.RemoteDataBaseBuilder;

/**
 * Allows the running of a discrete set of utilities in the context of the webserver,
 * providing various batch processing operations
 */
// To remove util functionality, comment out the line below
@WebServlet("/util/*")
public class UtilLauncher extends HttpServlet {
	private static final long serialVersionUID = 1L;
       
    /**
     * @see HttpServlet#HttpServlet()
     */
    public UtilLauncher() {
        super();
    }

	/**
	 * @see HttpServlet#doGet(HttpServletRequest request, HttpServletResponse response)
	 */
	protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException 
	{
		String util = request.getPathInfo().substring(1);
		System.out.println(util);

		OutputStreamWriter os = new OutputStreamWriter(response.getOutputStream());
		os.write("Running: "+util);
		switch (util)
		{
		case "DBB":
			try {
				DataBaseBuilder dbb = new DataBaseBuilder();
				dbb.run(os);
			} catch (SQLException e) {
				e.printStackTrace();
			}
			break;
		case "RDBB":
			try {
				RemoteDataBaseBuilder rdbb = new RemoteDataBaseBuilder();
				rdbb.run(os);
			} catch (SQLException e) {
				e.printStackTrace();
			}
			break;
		}
	    os.close();
	}

	/**
	 * @see HttpServlet#doPost(HttpServletRequest request, HttpServletResponse response)
	 */
	protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
	}

}
