package com.backend.collab;

import javax.servlet.annotation.WebServlet;
import org.eclipse.jetty.websocket.servlet.WebSocketServlet;
import org.eclipse.jetty.websocket.servlet.WebSocketServletFactory;
 
@SuppressWarnings("serial")
@WebServlet(name = "SyncServer", urlPatterns = { "/connect" })
public class SocketServlet extends WebSocketServlet
{
    @Override
    public void configure(WebSocketServletFactory factory) 
    {
        factory.getPolicy().setIdleTimeout(0);
        factory.register(WebSocketHandler.class);
    }

}