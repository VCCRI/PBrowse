package com.tests.logger;

import java.io.File;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

import javax.servlet.ServletException;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

/**
 * Servlet implementation class ResultsLogger
 */
@WebServlet("/testlogger")
public class ResultsLogger extends HttpServlet {
	private static final long serialVersionUID = 1L;
       
    /**
     * @see HttpServlet#HttpServlet()
     */
    public ResultsLogger() {
        super();
        // TODO Auto-generated constructor stub
    }

	/**
	 * @see HttpServlet#doGet(HttpServletRequest request, HttpServletResponse response)
	 */
	protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {}

	/**
	 * @see HttpServlet#doPost(HttpServletRequest request, HttpServletResponse response)
	 */
	protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException 
	{
		OutputStreamWriter os = new OutputStreamWriter(response.getOutputStream());
		String result = request.getParameter("result")+"\n\n";

		//create file if non existent
		new File("/opt/pbrowse/testresultlog.txt").createNewFile();
		
		try {
		    Files.write(Paths.get("/opt/pbrowse/testresultlog.txt"), result.getBytes(), StandardOpenOption.APPEND);
		}catch (IOException e) {
			os.write("There was an error, please contact the sysadmin.");
			os.close();
			
			e.printStackTrace();
			return;
		}
		
		os.write("Thankyou.");
		os.close();
	}

}
