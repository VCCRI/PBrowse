package com.backend.collab;

import java.util.logging.Level;
import java.util.logging.Logger;

import org.eclipse.jetty.websocket.api.Session;
import org.eclipse.jetty.websocket.api.annotations.OnWebSocketClose;
import org.eclipse.jetty.websocket.api.annotations.OnWebSocketConnect;
import org.eclipse.jetty.websocket.api.annotations.OnWebSocketError;
import org.eclipse.jetty.websocket.api.annotations.OnWebSocketMessage;
import org.eclipse.jetty.websocket.api.annotations.WebSocket;

import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.backend.collab.SessionManager;

import static com.backend.collab.SessionManager.getISM;

/**
 * Implements the websocket for maintaining shared sessions, passes received messages
 * onto the SessionManager instance after minimal parsing
 * @author root
 *
 */

@WebSocket(maxTextMessageSize = 64 * 1024)
public class WebSocketHandler 
{

	static Logger logger = Logger.getLogger(WebSocketHandler.class.getName());
	private final JsonParser jp = new JsonParser();
	
    @OnWebSocketConnect
    public void open(Session session) 
    {
    	getISM().openSocket(session);
    }
	
	@OnWebSocketClose
    public void onClose(Session session, int status, String reason) {
        System.out.println("Close: statusCode=" + status + ", reason=" + reason);
        getISM().closeSocket(session);
    }
	
	@OnWebSocketError
    public void onError(Throwable error) 
	{
		System.out.println("~~~~~~~~~~~~~~~~~~! SERIOUS ERROR !~~~~~~~~~~~~~~~~~~");
		error.printStackTrace();
	}
	
	/**
	 * The function invoked by the websocket whenever a new message is received
	 * @param message - the JSON encoded message
	 * @param session - the websocket session
	 */
	@OnWebSocketMessage
    public void handleMessage(Session session, String message) 
	{
		try
		{
			//parse the Stringified message into a Json object again
			JsonObject msg = jp.parse(message).getAsJsonObject();
            JsonObject params = msg.getAsJsonObject("params");
            
			//no message type specified, do nothing
            JsonElement ele = msg.get("type");
			if (ele == null) return;
			
			//all messages have a type
			String type = ele.getAsString();
			
			if (SessionManager.DEBUG)
			{
				logger.log(Level.INFO, "Request type: "+type+" "+msg.toString());
			}

			//filter TEST messages
			if (type.startsWith("TEST"))
			{
				getISM().testSync(session,type,params);
				return;
			}
			
			switch (type)
			{
			case "schat":
			case "update-status":
			case "update-tiers":
				if (params == null)
					return;
				
				getISM().notifySessionUsers(session, params);
				break;
				
			case "update-saved-config":
				if (params == null)
					return;
				
				getISM().updateSavedConfig(session,params);
				break;
				
			case "session-nominate-leader":
				if (params == null)
					return;
				
				getISM().setSessionLeader(session, params);
				break;
				
			case "register-remote-source":
				if (params == null)
					return;
				
				getISM().registerRemoteSource(session, params);
				break;
				
			case "session-blacklist-update":
				if (params == null)
					return;
				
				getISM().setSessionBlacklist(session, params);
				break;
				
			case "session-invite-user":
				if (params == null)
					return;
				
				getISM().doSessionInviteUser(session, params);
				break;
				
			case "session-kick-user":
				if (params == null)
					return;
				
				getISM().doSessionKickUser(session, params);
				break;
				
			case "reset-user-password":
				if (params == null)
					return;
				
				getISM().resetUserPassword(session, params);
				break;
				
			case "change-user-password":
				if (params == null)
					return;
				
				getISM().changeUserPassword(session, params);
				break;
				
			case "change-user-email":
				if (params == null)
					return;

				getISM().changeUserEmail(session, params);
				break;
				
			case "change-user-nickname":
				if (params == null)
					return;

				getISM().changeUserNickname(session, params);
				break;

			case "send-group-message":
				if (params == null)
					return;

				getISM().sendGroupEmail(session, params);
				break;

			case "leave-group":
				if (params == null)
					return;

				getISM().removeGroupUser(session, params);
				break;
				
			case "remove-group-file-acl":
				if (params == null)
					return;
				
				getISM().removeGroupFileACL(session, params);
				break;
				
			case "set-group-user-acl":
				if (params == null)
					return;
				
				getISM().updateGroupUserACL(session, params);
				break;
				
			case "remove-group-user":
				if (params == null)
					return;
				
				getISM().removeGroupUser(session, params);
				break;
				
			case "create-group":
				if (params == null)
					return;
				
				getISM().createGroup(session, params);
				break;
				
			case "delete-group":
				if (params == null)
					return;
				
				getISM().deleteGroup(session, params);
				break;
				
			case "add-group-user":
				if (params == null)
					return;
				
				getISM().addUserToGroup(session, params);
				break;
				
			case "share-group-file":
				if (params == null)
					return;
				
				getISM().shareGroupFile(session, params);
				break;
				
			case "get-group-info":
				if (params == null)
					return;
				
				getISM().getGroupInfo(session, params);
				break;
				
			case "get-user-groups":
				getISM().getUserGroups(session);
				break;
				
			case "get-all-comments":
				if (params == null)
					return;
				
				getISM().getAllComments(session, params);
				break;
				
			case "delete-comment":
				if (params == null)
					return;
				
				getISM().deleteComment(session, params);
				break;
				
			case "make-comment":
				if (params == null)
					return;
				
				getISM().makeComment(session, params);
				break;
				
            case "leader-options":
				if (params == null)
					return;
				
				getISM().sessionSetOption(session, params);
				break;
				
            case "resume-session":
            	if (params == null)
					return;
				
				getISM().doResumeSession(session, params);
				break;
				
			case "create-session":
				if (params == null)
					return;
				
				getISM().newSession(session,params);
				break;
				
			case "end-session":
				getISM().closeSession(session);
				break;
				
			case "join-session":
				if (params == null)
					return;
				
				getISM().joinSession(session, params);
				break;
				
			case "register-user":
				if (params == null)
					return;
				
				getISM().doRegister(session, params);
				break;
				
			case "get-user-files":
				getISM().getUserFiles(session);
				break;
				
			case "get-public-files":
				getISM().getPublicFiles(session);
				break;
				
			case "get-public-sessions":
				getISM().getPublicSessions(session);
				break;
				
			case "delete-user-file":
				if (params == null)
					return;
				
				getISM().deleteUserFile(session, params);
				break;
				
			case "toggle-user-file-public":
				if (params == null)
					return;
				
				getISM().setUserFilePublicStatus(session, params);
				break;
				
			case "login-user":
				if (params == null)
					return;
				
				getISM().doLogin(session, params);
				break;
				
			case "logout-user":
				getISM().doLogout(session);
				break;
				
			case "leave-session":
				getISM().leaveSession(session);
				break;
				
			default:
				break;
			}
		}
		catch (Exception e)
		{
			//logs the exceptions in /opt/tomcat8/logs/catalina.out
			logger.log(Level.SEVERE, null, e);
		}
		
	}

}