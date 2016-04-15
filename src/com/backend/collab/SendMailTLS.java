package com.backend.collab;

import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.Properties;

import javax.mail.Message;
import javax.mail.MessagingException;
import javax.mail.PasswordAuthentication;
import javax.mail.Session;
import javax.mail.Transport;
import javax.mail.internet.InternetAddress;
import javax.mail.internet.MimeMessage;

/**
 * Mail transmission module, maintains a mailing queue and sends messages asynchronously from
 * the main thread to prevent hang-ups. Currently uses a specially set-up gmail account to
 * forward messages.
 * @author root
 *
 */
public class SendMailTLS implements Runnable {

	private Session session;
	public boolean RUNNING = true;
	
	private ArrayList<MimeMessage> msgq;
	
	public SendMailTLS()
	{
		//read email configuration
		Properties p = new Properties();
		try {
			p.load(new FileInputStream("/opt/pbrowse/email_bot_config.properties"));
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		final String username = p.getProperty("email");
		final String password = p.getProperty("password");

		Properties props = new Properties();
		props.put("mail.smtp.auth", p.getProperty("mail.smtp.auth"));
		props.put("mail.smtp.starttls.enable", p.getProperty("mail.smtp.starttls.enable"));
		props.put("mail.smtp.host", p.getProperty("mail.smtp.host"));
		props.put("mail.smtp.port", p.getProperty("mail.smtp.port"));

		//initialize the mail instance
		session = Session.getInstance(props,
			new javax.mail.Authenticator() {
			protected PasswordAuthentication getPasswordAuthentication() {
			return new PasswordAuthentication(username, password);
			}
		});
		
		msgq = new ArrayList<MimeMessage>();
	}
	
	/**
	 * Add a new message to the sending queue
	 * @param subject - email subject
	 * @param to - to address
	 * @param msg - message html
	 */
	public void enqueueMessage(String subject, String to, String msg) 
	{
		try {
			MimeMessage message = new MimeMessage(session);
			message.setFrom(new InternetAddress("noreply@pbrowse.com"));
			message.setRecipients(Message.RecipientType.TO,InternetAddress.parse(to));
			message.setSubject(subject);
			message.setText(msg, "utf-8", "html");
			
			msgq.add(message);
		} catch (MessagingException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * The mailing thread, checks for new messages to send every 3 seconds and sends them
	 */
	@Override
	public void run() 
	{
		while (this.RUNNING)
		{
			try {
				Thread.sleep(3000);
			
				if (!msgq.isEmpty())
				{
					MimeMessage message = msgq.remove(0);
					Transport.send(message);
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
}
