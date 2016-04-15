/**
 * The CollabManager is the super-object behind all communications and synchronisation features
 * of PBrowse. It works by passing predetermined messages to the host server, with customisable
 * parameters. The manager waits for a response from the server and processes it accordingly.
 * Only 1 instance per user is created.
 */
function CollabManager () 
{
	this.sendMessage = function(m)
	{
//		console.info(m);
		this.sock.send(m);
	};
	
	//datatables
	this.dd_private_data = {};
	this.dd_public_data  = {};
	this.dd_group_data 	 = {};
	
	//privacy constants
	this.PRIVACY_PUBLIC = 4;
	this.PRIVACY_PRIVATE = 2;
	this.PRIVACY_ONLYME = 1;
	
	//id to type
	this.id_to_type = ["bam","bigbed","bigwig","vcf","bed","wig","2bit","das"];
	
	//pbrowse saved tier configurations
	this.saved_config = {};
	
	/**
	 * Storage of variables used by PBrowse
	 */
	this.public_sessions_list = {};
	this.session_blacklist = {};
	this.active_session = undefined;
	this.active_genome = "hg38"; //default to human
	this.me = undefined;
	this.grouplist = {};
	this.uauth = undefined;
	this.session_users = {};
	this.comments = {};
	this.tierIDMap = {};
	this.group_mem_list = {};
	this.group_shared_files = {};
	this.active_group = undefined;
	this.rootComment = null;
	this.update_timeout = null;
	this.tiers_timeout = null;
	this._follower_diverge = true;
	this._follower_sync = true;
	this._private_session = true;
	this._numloadedtracks = 0;
	this.dialog_stack = [];
	
	/**
	 * Adds the file specified by the configuration 'json' and identified by the ID 'id'
	 * If and only if the file is reachable by the current user.
	 */
	this._addTrack = function(id, json)
	{
		var http = new XMLHttpRequest();
	    http.open('HEAD', "//"+window.location.host+"/files/"+id);
	    http.onreadystatechange = function() 
	    {
	        if (this.readyState == this.DONE) 
	        {
	            if (this.status == 200 || this.status == 206)
            	{
	            	b.addTier(json);
            	}
	            else
            	{
	            	//console.error("User lacks permissions to view selected file. ERR: "+this.status);
            	}
	        }
	    };
	    http.send();
	};
	/**
	 * Removes the comment 'c' from the comment tree, trimming the parent node as well, if
	 * it has no more children
	 */
	this._removeComment = function(c)
	{
		if (c == undefined)
			return;
		
		if (c.isnew)
		{
			c.setnew(-1);
		}
		c.deleteme();
		
		delete cm.comments[c.id];

		var parent = c.parent;
		delete parent.children[c.id];
		if (Object.keys(parent.children).length == 0)
		{
			parent.deleteme();
			
			//update trackref node
			var tnode = parent.parent;
			delete tnode.children[parent.chr];
		}
	};
	/**
	 * Returns TRUE if the current user is the leader of their collaborative session
	 */
	this._amILeader = function()
	{
		var me = this.session_users[this.getUser()];
		if (me == undefined || me == null)
		{
			return false;
		}
		return me.isleader;
	};
	/**
	 * Return true if the user has currently got the file with id loaded
	 */
	this._viewingFile = function(id)
	{
		var r = -1;
		var tiers = this.getTiers();
		for (var t=0; t<tiers.length; t++)
		{
			if (tiers[t].source["_id"] == id)
			{
				r = t;
			}
		}
		return r;
	};
	/**
	 * Toggles a mask on the followers browser preventing further interaction with it.
	 */
	this._allowDivergance = function(s)
	{
		if (!this._amILeader())
		{
			return;
		}
		
		this._follower_diverge = s;

		var p = {
			uuid: this.active_session,
			user: this.getUser(),
			setoption: "allowdivergence:"+s+";"
		};

		var message = {
			type: "leader-options",
			params: p
		};
		
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Sets the session as public or private depending on the parameter 's'
	 */
	this._setPrivateSession = function(s)
	{
		if (!this._amILeader())
		{
			return;
		}
		
		this._private_session = s;

		var p = {
			uuid: this.active_session,
			user: this.getUser(),
			setoption: "private:"+s+";"
		};

		var message = {
			type: "leader-options",
			params: p
		};
		
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Toggles the ability of all followers to diverge views from the leaders view
	 * regardless of whether the leader is still actively changing it
	 */
	this._enableFollowerSync = function(s)
	{
		if (!this._amILeader())
		{
			return;
		}
		
		this._follower_sync = s;

		var p = {
			uuid: this.active_session,
			user: this.getUser(),
			setoption: "followersync:"+s+";"
		};

		var message = {
			type: "leader-options",
			params: p
		};
		
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Removes any tiers the calling user does not directly own, making sure that if they
	 * log out, dalliance will not attempt to load them and fail.
	 */
	this._removePrivateTiers = function()
	{
		var t = this.getTiers();
		for (var i=0; i<t.length; i++)
		{
			//remove if private tier and the user isn't the owner of the file
			if (t[i]["source"]["_ispublic"] != this.PRIVACY_PUBLIC) 
			{
				if (t[i]["source"]["_owner"] != undefined)
				{
					if (t[i]["source"]["_owner"] != this.getUser())
					{
						b.removeTier({"index": i});
					}
				}
			}
		}
	};
	/**
	 * Returns the USERNAME of the calling user, as stored in cookies
	 * The user must be logged in.
	 */
	this.getUser = function ()
	{
		if (this.me != undefined)
			return this.me;
		
		var val = undefined;
		var cs = document.cookie.split(";");
		for (var c=0; c<cs.length; c++)
		{
			if (cs[c].indexOf("username") != -1)
			{
				val = cs[c].split("=")[1];
				break;
			}
		}
		this.me = val;
		
		return val;
	};
	/**
	 * Returns the UAUTH token of the current user used to verify the caller's identity
	 * User must be logged in.
	 */
	this.getUAuth = function ()
	{
		if (this.uauth != undefined)
			return this.uauth;
		
		var val = undefined;
		var cs = document.cookie.split(";");
		for (var c=0; c<cs.length; c++)
		{
			if (cs[c].indexOf("uauth") != -1)
			{
				val = cs[c].split("=")[1];
				break;
			}
		}
		this.uauth = val;
		
		return val;
	};
	/**
	 * Returns the current set of tiers loaded by the calling user's genome browser, along with
	 * their configuration objects.
	 */
	this.getTiers = function ()
	{
		var conf = [];
		for (var i=0; i<b.tiers.length; i++)
		{
			conf.push({"source": b.tiers[i].dasSource, "config": b.tiers[i].config || {}});
		}
		return conf;
	};
	/**
	 * Loads the tiers configured by the 'conf' parameter in place of any currently loaded tiers
	 */
	this.setTiers = function (conf)
	{
		var remapping = [];
		//prefill with invalid entries, if the track is not found in the
		//configuration conf, it will be removed
		for (var i=0; i<b.tiers.length; i++)
		{
			remapping[i] = -1;
		}
		
		//check all defined config styles are the same, otherwise the
		//tier must be re-styled
		var stylesEqual = function(a,b)
		{
			var k = Object.keys(a);
			var l = Object.keys(b);
			
			if (k.length != l.length)
				return false;
			
			for (var s=0; s<k.length; s++)
			{
				if (a[k[s]] != b[k[s]])
					return false;
			}
			return true;
		};
		
		//array of files for which comments must be loaded
		var newfiles = [];
		
		for (var i=0; i<conf.length; i++)
		{
			//check if we have already loaded comments for this file
			//if so, we do not waste time loading them again
			var id = conf[i].source["_id"];
			if (id != undefined)
			{
				if (this.rootComment.children[id] == undefined)
				{
					newfiles.push(id);
				}
			}
			
			//search for the tier matching the config object to see if it is
			//already loaded
			var sourceFound = false;
			for (var j=0; j<b.tiers.length; j++)
			{
				//sources are equivalent - the tracks are the same
				if (sourcesAreEqual(b.tiers[j].dasSource,conf[i].source) )
				{
					//if the styles have been changed, update our copy
					if (!stylesEqual(b.tiers[j].config,conf[i].config))
					{
						b.tiers[j].config = conf[i].config;
						b.tiers[j]._updateFromConfig();
					}
					
					//create mapping of position for reordering tracks
					//move index i -> map[i]
					remapping[j] = i;
					sourceFound = true;
					break;
				}
			}
			if (!sourceFound)
			{
				//source not found in existing set? add it
				if (conf[i].source["_ispublic"] == this.PRIVACY_ONLYME)
				{
					if (conf[i].source["_owner"] != this.getUser())
					{
						//don't render "onlyme" tracks if not owner - in any circumstance
						continue;
					}
				}
				
				//set new remapping to put the tier in correct place
				remapping[b.tiers.length] = i;
				
				//make the tier from config information
				var tier = b.makeTier(conf[i].source,conf[i].config || {});
				b.refreshTier(tier);
			}
		}

		//do tier arrangement
		var tset = [];
		for (var i=0; i<b.tiers.length; i++)
		{
			tset.push(b.tiers[i]);
		}
		for (var j=0; j<remapping.length; j++)
		{
			//shuffle tier
			if (remapping[j] >= 0)
				b.tiers[remapping[j]] = tset[j];
			//delete tier
			else
			{
				b.nn_removeTier({"index": j});
			}
		}
		//draw result
		b.reorderTiers();
		
		//fix refresh removing everything stupidly
		b.storeTierStatus();
		
		//load comments only if we haven't loaded them once
		this.getViewedFiles();
		this.cleanComments();
		this.getAllComments(newfiles);
	};
	/**
	 * Returns the position status of the calling user and any highlights they have made
	 */
	this.getStatus = function ()
	{
		var username = this.getUser();
		if (username === undefined)
		{
			return undefined;
		}
		var status = {
			"username" : username,
			"chr" : b.chr,
			"start_win" : b.viewStart|0,
			"end_win" : b.viewEnd|0,
			"highlights" : b.highlights,
			"genome" : this.active_genome
		};
		return status;
	};
	/**
	 * Attempts to join the session identified by the 'id' parameter with the password
	 * provided by the 'code' parameter. Sends the 'join-session' message to the server.
	 */
	this.joinSession = function (id,code) 
	{
		var status = this.getStatus();
		if (status === undefined)
		{
			err($("#joinsessionError"),"Only logged in users may join a collaborative session.","danger");
			return;
		}
		if (this.active_session)
		{
			err($("#joinsessionError"),"You must leave your active session before joining another.","danger");
			return;
		}
		if (id === null || id === undefined || id === "")
		{
			err($("#joinsessionError"),"You must enter the id of the session to join.","danger");
			return;
		}
		
		status["code"] = code;

		status.uuid = id;
		var message = {
			type: "join-session",
			params: status
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Creates a new session if the user is not already part of a session. Sends the
	 * 'create-session' message to the server.
	 */
	this.newSession = function ()
	{
		if (this.getUser() === undefined)
		{
			showAlert("Error","Only logged in users may start a collaborative session.");
			return;
		}
		if (this.active_session)
		{
			showAlert("Error","You already belong to a session.");
			return;
		}
		
		var isprivate = document.getElementById("newSessionPrivate").checked;
		document.getElementById("checkPrivateSession").checked = isprivate;
		
		var code = $("#newSessionCode").val(),
			name = $("#newSessionName").val();
		
		var p = {
			"private" : isprivate,
			"code" : code,
			"name" : name,
			"status" : this.getStatus(),
			"tiers" : this.getTiers(),
			"genome" : this.active_genome,
		};
		
		if (this.active_genome == "generic")
		{
			p["custom_genome"] = genomes["generic"].coordSystem.speciesName;
		}

		//set session label and code text if provided by the user
		$("#psessionCode").text($("#newSessionCode").val()!=""?" | Code: "+code:"");
		$("#sessionNameLabel").text(name);
		
		var message = {
			type: "create-session",
			params: p,
		};
		this.sendMessage(JSON.stringify(message));

//		this._updateStatus();
	};
	/**
	 * Sends the user's current status to all participating session users - invoked by 
	 * user actions only
	 */
	this.updateStatus = function()
	{
		this.update_timeout = null;
		
		var status = this.getStatus();
		if (status != null)
		{
			status.human = true;
			var message = {
				type: "update-status",
				params: status
			};
			this.sendMessage(JSON.stringify(message));
		}
	}.bind(this);
	/**
	 * The machine invoked update request - same as previous, but it will not
	 * cause further update propagation
	 */
	this._updateStatus = function()
	{
		var status = this.getStatus();
		if (status != null)
		{
			status.human = false;
			var message = {
				type: "update-status",
				params: status
			};
			this.sendMessage(JSON.stringify(message));
		}
	}.bind(this);
	/**
	 * Sends a session chat message to all participating session members
	 */
	this.sendSessionMsg = function(msg, raw)
	{
		if (raw == undefined)
			raw = false;
		
		var p = {
			"smsg" : msg,
			"raw" : raw,
		};

		var message = {
			type: "schat",
			params: p,
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Removes the calling user from their session. If they are the leader, the role is
	 * automatically passed on to another user. Has no effect if the user is not part of
	 * a session
	 */
	this.leaveSession = function()
	{
		//notification text
		this.sendSessionMsg(this.getUser()+" has left the session.", true);
		
		var message = {
			type: "leave-session"
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Retrieves a list of the currently available public collaborative sessions
	 */
	this.getPublicSessions = function()
	{
		var message = {
			type: "get-public-sessions"
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Retrieves a list of all the files the calling user has uploaded to PBrowse
	 */
	this.getUserFiles = function()
	{
		if (this.getUser() == undefined)
		{
			return;
		}
		
		var message = {
			type: "get-user-files"
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Retrieves a list of all the files uploaded by all the users of PBrowse who have specified
	 * them as publically accessible.
	 */
	this.getPublicFiles = function()
	{
		var message = {
			type: "get-public-files"
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Allows the leader of a collaborative session to instantly end a collaborative session and
	 * remove all of its users.
	 */
	this.endSession = function()
	{
		var message = {
			type: "end-session",
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Registers a new user account using the provided parameters. 
	 */
	this.doRegister = function(username, nickname, email, password)
	{
        var p = {
            "username" : username,
            "nick" : nickname,
            "email" : email,
            "password" : password
        };

		var message = {
			type: "register-user",
            params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Attempts to login the the user identified by the 'username' and 'password' parameter
	 */
	this.doLogin = function(username, password)
	{
        var p = {
            "username" : username,
            "password" : password
        };

		var message = {
			type: "login-user",
            params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Attempts to delete the file identified by the 'id' parameter. If the invoker does not
	 * own the file it will fail and no change is made.
	 */
	this.deleteDataFile = function(id)
	{
        var p = {
            "id" : id,
        };

		var message = {
			type: "delete-user-file",
            params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Swithes the public status of the file identified by the 'id' parameter from public to
	 * private, or vice versa.
	 */
	this.togglePublicDataFile = function(id, status)
	{
        var p = {
            "id" : id,
            "privacy" : status,
        };

		var message = {
			type: "toggle-user-file-public",
            params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Sends an update to all the members of collaborative session informing them of a change
	 * to the calling user's track configuration.
	 */
	this._updateTiers = function()
	{
		var p = {
            "tiers" : this.getTiers(),
            "genome" : this.active_genome,
        };
		
		var custom_genome = undefined;
		if (this.active_genome == "generic")
		{
			custom_genome = genomes["generic"].coordSystem.speciesName;
		}
		p["custom_genome"] = custom_genome;

		var message = {
			type: "update-tiers",
            params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	this.updateTiers = this._updateTiers.bind(this);
	/**
	 * Attempts to log out the calling user. Has no effect if the user is not logged in.
	 */
	this.doLogout = function()
	{
		var message = {
			type: "logout-user",
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Processes the 'opts' parameter for predetermined arguments and applies any changes
	 * to the configuration of the user's collaborative session. Only the session leader can
	 * invoke changes, followers must accept any configuration changes they receive.
	 */
	this.processOptions = function(opts)
	{
		this._follower_sync = opts["followersync"]=="true"?true:false;
		
		//update option box for when user ever becomes leader
		document.getElementById("checkFollowerSync").checked = this._follower_sync;
		
		//show session code and name
		document.getElementById("psessionCode").innerHTML = "";
		document.getElementById("psessionCode").appendChild(
				document.createTextNode(opts["code"]!=""?" | Code: "+opts["code"]:"")
				);
		document.getElementById("sessionNameLabel").innerHTML = "";
		document.getElementById("sessionNameLabel").appendChild(
				document.createTextNode(opts["name"])
				);
		
		this._private_session = opts["private"]=="true"?true:false;
		document.getElementById("checkPrivateSession").checked = this._private_session;
		
		if (this._amILeader() == false)
		{
			this._follower_diverge = opts["allowdivergence"]=="true"?true:false;
			document.getElementById("checkFollowerDiverge").checked = this._follower_diverge;
			
			//create an input capturing, transparent mask over the biodalliance browser
			if (!this._follower_diverge)
			{
				$("#followerDivergeMask").show();
				$("#followerDivergeMask").focus();
			}
			else $("#followerDivergeMask").hide();
		}
		else
		{ 
			this._follower_diverge = true;
			$("#followerDivergeMask").hide();
		}
	};
	/**
	 * Makes a new comment on the file identified by the 'trackref' parameterm using the other
	 * provided parameters. Public comments are sent to every active user viewing the same file,
	 * regardless of whether they are in a collaborative session together or not.
	 */
	this.makeComment = function(trackref, chr, start, end, ctext, ispublic)
	{
        var p = {
            "trackref" : trackref,
            "chr" : chr,
            "start" : start,
            "end" : end,
            "ctext": ctext,
            "ispublic": ispublic
        };

		var message = {
			type: "make-comment",
            params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	this.cleanComments = function()
	{
		//trim comments for files we don't see
		var rem = [];
		for (var i in this.comments)
		{
			if (this.tierIDMap[this.comments[i].trackref] == undefined)
			{
				rem.push(i);
			}
		}
		
		//delete all the comments we shouldn't see
		for (var i=0; i<rem.length; i++)
		{
			this._removeComment(this.comments[rem[i]]);
		}
		
		if (this.rootComment !== null)
		{
			//for the comment tree simply delete the trackref level nodes
			rem = [];
			for (var c in this.rootComment.children)
			{
				if (this.tierIDMap[c] == undefined)
				{
					rem.push(c);
				}
			}
			
			//delete all the trackref nodes we shouldn't see
			for (var i=0; i<rem.length; i++)
			{
				//remove from dom
				this.rootComment.children[rem[i]].deleteme();
				delete this.rootComment.children[rem[i]];
			}
		}
	};
	this.getViewedFiles = function()
	{
		this.tierIDMap = {};

		var t = this.getTiers();
		for (var i=0; i<t.length; i++)
		{
			if (t[i]["source"]["_id"] !== undefined)
			{
				//record mapping of ids to tier names for human readibility
				this.tierIDMap[t[i]["source"]["_id"]] = t[i]["source"]["name"];
			}
		}
	};
	
	/**
	 * Retrieves a list of all the comments made on the files loaded as tracks in the browser, 
	 * which have been made by the calling user - or made publically viewable. Removes comments
	 * for files which the user is no longer viewing
	 */
	this.loadComments = function()
	{
		//determine which tiers we must load comments for
		this.getViewedFiles();
		this.cleanComments();
		
		//do load comments for all appropriate tracks
		var arr = Object.keys(this.tierIDMap);
		if (arr.length > 0)
			this.getAllComments(arr);
	};
	/**
	 * Loads all comments for the files specified in the 'trackref' array parameter
	 */
	this.getAllComments = function(trackref)
	{
        var p = {
            "trackref" : trackref,
        };

		var message = {
			type: "get-all-comments",
            params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Deletes the comment identified by the 'id' and 'trackref' parameters. Only the 
	 * creator of a comment has the right to delete it.
	 */
	this.deleteComment = function(c)
	{
		var p = {
            "id" : c.id,
            "trackref" : c.trackref,
        };

		var message = {
			type: "delete-comment",
            params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Creates a new user group using the provided 'name' parameter. It will fail if the
	 * name has already been taken
	 */
	this.createGroup = function(name)
	{
		var p = {
            "groupname" : name,
        };

		var message = {
			type: "create-group",
            params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Retrieves a list of all the groups to which the calling user belongs
	 */
	this.getUserGroups = function()
	{
		var message = {
			type: "get-user-groups",
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Retrieves a list of all the users belonging to the group identified by the 'groupname'
	 * parameter. A group must have at least 1 member - the creator.
	 */
	this.getGroupUsers = function(groupname)
	{
		var p = {
            "groupname" : groupname,
        };

		var message = {
			type: "get-group-users",
            params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Retrieves a list of all the files shared with the group identified by the 'groupname'
	 * parameter, as well as all the users in that group.
	 */
	this.getGroupInfo = function(groupname)
	{
		var p = {
            "groupname" : groupname,
        };

		var message = {
			type: "get-group-info",
            params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Attmpts to share a file identified by the 'fileid' parameter to the group identified by
	 * the 'groupname' parameter. It will fail if the calling user lacks the appropriate sharing
	 * permission on the group.
	 */
	this.shareFileToGroup = function(fileid, groupname)
	{
		var p = {
			"fileid" : fileid,
            "groupname" : groupname,
        };

		var message = {
			type: "share-group-file",
            params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Attempts to add the user identified by the 'username' parameter to the group identified
	 * by the 'groupname' parameter, with an intial set of priileges according to the 'access_level'
	 * parameter. It will fail if the calling user lacks the privilege to manage group users.
	 */
	this.addUserToGroup = function(username, groupname, access_level)
	{
		var p = {
			"username" : username,
            "groupname" : groupname,
            "access_level" : access_level,
        };

		var message = {
			type: "add-group-user",
            params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Attempts to delete the group idenitifed by the 'groupname' parameter. It fails if the calling
	 * user is not the owner of the group.
	 */
	this.deleteGroup = function(groupname)
	{
		var p = {
	        "groupname" : groupname,
	    };
	
		var message = {
			type: "delete-group",
	        params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Attempts to modifiy the permission set of the identified user, in the identified group,
	 * according to the permissions set in the 'newacl' parameter. Fails if the calling user tries
	 * to modify their own permissions, or grant permissions they themselves have no privilege to
	 * grant
	 */
	this.updateGroupUserACL = function(groupname, username, newacl)
	{
		var p = {
	        "groupname" : groupname,
	        "username" : username,
	        "acl" : newacl
	    };
	
		var message = {
			type: "set-group-user-acl",
	        params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Removes the user from the group if the calling user has the appropriate group user
	 * management permission.
	 */
	this.removeUserFromGroup = function(groupname, username)
	{
		var p = {
	        "groupname" : groupname,
	        "username" : username,
	    };
	
		var message = {
			type: "remove-group-user",
	        params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Updates the access requirement of the given file, in the specified group with the new
	 * access level defined by the 'newacl' parameter. NOW DEPRECATED.
	 */
	this.updateGroupFileACL = function(groupname, fileid, newacl)
	{
		var p = {
	        "groupname" : groupname,
	        "fileid" : fileid,
	        "acl" : newacl
	    };
	
		var message = {
			type: "update-group-file-acl",
	        params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Un-shares the file with the group, preventing members of that group from accessing it
	 */
	this.removeGroupFileACL = function(groupname, fileid)
	{
		var p = {
	        "groupname" : groupname,
	        "fileid" : fileid
	    };
	
		var message = {
			type: "remove-group-file-acl",
	        params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Removes the calling user from the specified group. Given valid parameters, the call will
	 * never fail.
	 */
	this.leaveGroup = function(groupname)
	{
		var p = {
	        "groupname" : groupname
	    };
	
		var message = {
			type: "leave-group",
	        params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Sends an email message to all the members of the group 'groupname' with message text
	 * 'message'. Only the owner can invoke this function successfully. Otherwise it will silently
	 * fail.
	 */
	this.sendGroupMesssage = function(groupname, message)
	{
		var p = {
	        "groupname" : groupname,
	        "message" : message
	    };
	
		var message = {
			type: "send-group-message",
	        params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Invokes a password change for the calling user, succeeds only if the old password is correct
	 * and the new passwords match.
	 */
	this.changeUserPassword = function(old, new1, new2)
	{
		if (old.length == 0 || new1.length == 0 || new2.length == 0)
			return;
		
		var p = {
	        "old" : old,
	        "password1" : new1,
	        "password2" : new2,
	    };
	
		var message = {
			type: "change-user-password",
	        params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Changes the calling user's nickname. Will always succeed.
	 */
	this.changeUserNickname = function(nick)
	{
		if (nick.length == 0)
			return;
		
		var p = {
	        "nickname" : nick,
	    };
	
		var message = {
			type: "change-user-nickname",
	        params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Will verify the email address by generating a token and sending a message to the new
	 * address. If it is retrieved by the user, then the email is updated.
	 */
	this.changeUserEmail = function(email)
	{
		var p = {
	        "email" : email,
	    };
	
		var message = {
			type: "change-user-email",
	        params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Initiates a password reset request, using either the username or email of an existing user
	 */
	this.doPasswordReset = function()
	{
		var p = {
	        "usermail" : $("#l_username").val(),
	    };
	
		var message = {
			type: "reset-user-password",
	        params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Sets the provided user as leader - provided the caller is leader
	 */
	this.sessionNominateLeader = function(leader)
	{
		var p = {
	        "leader" : leader,
	    };
	
		var message = {
			type: "session-nominate-leader",
	        params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Prevents a target user from joining the session - provided the caller is leader
	 */
	this.sessionBlacklistUpdate = function(user, option)
	{
		var p = {
	        "user" : user,
	        "option" : option,
	    };
	
		var message = {
			type: "session-blacklist-update",
	        params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Invites the targeted user to join this session - provided the caller is leader
	 */
	this.sessionInviteUser = function(user)
	{
		var p = {
	        "user" : user,
	    };
	
		var message = {
			type: "session-invite-user",
	        params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Forcibly removes the targetted user from the session - provided the caller is leader
	 */
	this.sessionKickUser = function(user)
	{
		var p = {
	        "user" : user,
	    };
	
		var message = {
			type: "session-kick-user",
	        params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	this.registerRemoteSource = function(
			dataurl, data_format, feature_seq, studyid,
			trackname, genometag, description, ispublic
			)
	{
		var p = {
	        "dataurl" : dataurl,
	        "data_format" : data_format,
	        "feature_seq" : feature_seq,
	        "studyid" : studyid,
	        "trackname" : trackname,
	        "genometag" : genometag,
	        "description" : description,
	        "ispublic" : ispublic
	    };
	
		var message = {
			type: "register-remote-source",
	        params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Draws the new blacklist as a tooltip for the glyph
	 */
	this.drawBlacklist = function()
	{
		//update blacklist tooltip
		var str = "";
		var names = Object.keys(this.session_blacklist);
		if (names.length > 0)
		{
			for (var e=0; e<names.length; e++)
			{
				str += names[e]+" ";
			}
			$("#blacklistTooltip").attr("data-original-title",str);
		}
		else
		{
			$("#blacklistTooltip").attr("data-original-title","<Empty>");
		}
	};
	/**
	 * Updates the stored configuration lists for the current user
	 */
	this.updateSavedConfig = function()
	{
		var cfg = {};
		for (var n in this.saved_config)
		{
			var scfg = this.saved_config[n]; 
			cfg[scfg.name] = 
			{
				"name": scfg.name, "desc":scfg.desc, 
				"genome":scfg.genome, "time":scfg.time, 
				"tiers":scfg.tiers
			};
		}
		
		var p = {
	        "config" : JSON.stringify(cfg),
	    };
	
		var message = {
			type: "update-saved-config",
	        params: p
		};
		this.sendMessage(JSON.stringify(message));
	};
	/**
	 * Updates current user's view to match the provided genome and tier set
	 */
	this.synchGenomes = function(genome,tiers,custom_genome)
	{
		var _t = this;
		
		//fix for synchronising custom genomes
		if (genome != undefined)
		{
			if (genome != this.active_genome)
			{
				if (genome == "generic")
				{
					var gen = genomes["generic"];
					gen.coordSystem.speciesName = custom_genome;
					gen.sources = [];
					for (var i=0; i<tiers.length; i++)
					{
						gen.sources.push(tiers[i].source);
					}
				}
				this.changeGenome(genome, function(){
					b.nn_removeAllTiers();
					_t.setTiers(tiers);
				});
			}
			else
			{
				_t.setTiers(tiers);
			}
		}
	};
	/**
	 * Function responsible for listening to and processing responses from the PBrowse server
	 * and implementing changes in the client
	 */
	this.onmessage = function (event)
	{
		var msg = JSON.parse(event.data);
		switch (msg.rtype)
		{
		
			/**
			 * A new session was successfully created by the calling user
			 */
			case "new-session":
				$('#newSessionDLG').modal('hide');
				
				//disable creating or joining more sessions while controlling a session
				$('#cs_menu').parent().addClass("disabled");
				$('#js_menu').parent().addClass("disabled");
				
				$("#sessionPanel").removeClass("hidden",250);
				
				this.active_session = msg.session;
				$("#activeSessionDiv").show();
				$("#adminControllerDiv").show();
				$("#sessionLeaderIndicator").text("| Leader");
				$("#sessionChatPanel").removeClass("hidden");
				$("#psessionID").text(msg.session);
				
				b.storeTierStatus();
				break;

			/**
			 * The calling user either closed, or left the session, for the client, the result
			 * is the same in both cases
			 */
			case "close-session":
			case "leave-session":
				if (msg.username === this.getUser() || this.getUser() === undefined || msg.rtype == "close-session")
				{
					//reenable create/join
					if (this.getUser() != undefined)
					{
						$('#cs_menu').parent().removeClass("disabled");
						$('#js_menu').parent().removeClass("disabled");
					}

					if (!this._amILeader() && msg.rtype == "close-session")
					{
						showAlert("Error","The session was terminated by the leader");
					}
					
					this.active_session = undefined;
					this.session_users = {};
					this._removePrivateTiers();
					$("#sessionPanel").addClass("hidden",250);
					
					$("#activeSessionDiv").hide();
					$("#adminControllerDiv").hide();
					$("#sessionChatPanel").addClass("hidden");
					
					$("#userListDiv").html("");
					
					//remove non-owned private files from file listing
					if (this.getUser() != undefined)
					{
						var _t = this;
						Object.keys(_t.dd_private_data).forEach(function (i){
							if (_t.dd_private_data[i].owner != _t.getUser())
							{
								_t.dd_private_data[i].deleteme();
								delete _t.dd_private_data[i];
							}
						});
					}

					this._follower_diverge = true;
					this._follower_sync = true;
					
					$("#followerDivergeMask").hide();
					
					$("#psessionCode").html("");
					b.storeTierStatus();
				}
				else
				{
					//another user has left the session - not ME, remove him from
					//our user lists
					this.session_users[msg.username].deleteme();
					delete this.session_users[msg.username];

					this.updateUserList(msg["user-status"]);
				}
				break;
				
			/**
			 * The user attempted to join a collaborative session
			 */
			case "join-session":
				if (msg.status == "connected")
				{
					$("#sessionPanel").removeClass("hidden",250);
					$("#activeSessionDiv").show();

					//process file list -- always process this
					for (var i=0; i<msg.allfiles.length; i++)
					{
						//make sure its not already in our table
						if (this.dd_private_data[msg.allfiles[i].id] == null)
						{
							//cannot see "onlyme" files
							if (msg.allfiles[i].ispublic == this.PRIVACY_ONLYME)
								continue;
							
							//this is another user's file which we will temporarily appropriate
							this.filetable.row.add([ msg.allfiles[i].id, msg.allfiles[i].owner, msg.allfiles[i].trackname,
							                         msg.allfiles[i].description, msg.allfiles[i].path, msg.allfiles[i].format, msg.allfiles[i].ispublic, 
							                         msg.allfiles[i].indexedby, msg.allfiles[i].isindex, msg.allfiles[i].remote])
							                         .draw( false );
						}
					}
					this.updateUserList(msg["user-status"]);

					//all session users are notified of this joining
					if (msg.joiner != this.getUser())
						return;
					
					$('#joinSessionDLG').modal('hide');
					
					this.active_session = msg.session;
					
					//disable creating or joining more sessions while in a session
					$('#cs_menu').parent().addClass("disabled");
					$('#js_menu').parent().addClass("disabled");
					
					//process session options
					this.processOptions(msg.options);
					
					//update tiers to leader's (if this is not a returning leader)
					//fix for synchronising custom genomes
					this.synchGenomes(msg.genome, msg.tiers, msg["custom_genome"]);
					
					//process chat history
					if (msg.chathistory != undefined)
					{
						var schat = document.getElementById('schatArea');
						for (var i=0; i<msg.chathistory.length; i++)
						{
							schat.appendChild(document.createTextNode(msg.chathistory[i]));
							schat.appendChild(document.createElement("br"));
						}
						var js = $(schat);
						if(js.length)
							js.scrollTop(js[0].scrollHeight - js.height());
					}
					
					//update the user's copy of the session user blacklist
					if (msg.blacklist != undefined)
					{
						this.session_blacklist = {};
						
						for (var i=0; i<msg.blacklist.length; i++)
						{
							this.session_blacklist[msg.blacklist[i]] = 0;
						}
						this.drawBlacklist();
					}
					else
					{
						this.session_blacklist = {};
					}
					
					//switch view back to genome viewer if required
					$("#public-session-dlg").modal("hide");
					
					$("#psessionID").html(this.active_session);
					$("#sessionChatPanel").removeClass("hidden");
					
					//notification text
					this.sendSessionMsg(this.getUser()+" has joined the session.", true);
				}
				else
				{
					err($("#joinsessionError"),msg.error,"danger");
				}
				break;
			
			/**
			 * The list of public sessions as provided by the server, the client assembles
			 * them into a table which the user can interact with
			 */
			case "public-sessions":
				this.publicsessiontable.clear().draw();
				for (var i=0; i<msg.sessions.length; i++)
				{ 
					//add the session entry as a row in the public session table
					this.publicsessiontable.row.add([msg.sessions[i].uuid,msg.sessions[i].name,
					                                 msg.sessions[i].numusers,msg.sessions[i].code])
					                                 .draw();
				}
				break;
				
			/**
			 * The user attempted to delete a file
			 */
			case "file-delete":
				if (msg.status == "success")
                {
					//if the file is listed in public/private file tables, remove it
					$( "#file-manager-dlg" ).modal("hide");
					if (this.dd_private_data[msg.id])
					{
						this.dd_private_data[msg.id].deleteme();
						delete this.dd_private_data[msg.id];
					}
					if (this.dd_public_data[msg.id])
					{
						this.dd_public_data[msg.id].deleteme();
						delete this.dd_public_data[msg.id];
					}
					
					//check if we are viewing the file - and remove it if necessary
					var id = this._viewingFile(msg.id);
					if (id != -1)
					{
						showAlert("File Permissions", "File: "+b.tiers[id].dasSource.name+" was deleted.");
						b.nn_removeTier({"index": id});
					}
                }
				else
				{
					err($("#fmdlg_response"),msg.error,"danger");
				}
				break;
				
			/**
			 * A new private file was uploaded by a user - either in a session or not
			 */
			case "new-session-file":
				if (msg.remote == true)
				{
					err($("#remoteuploadError"),"File registered with PBrowse!","success");
				}
				
				this.filetable.row.add([ msg.id, msg.owner, msg.trackname,
				                   		 msg.description, msg.path, msg.format, msg.ispublic, 
				                   		 msg.indexedby, msg.isindex, msg.remote])
				                         .draw( false );
				break;
				
			/**
			 * The leader sent a group email
			 */
			case "send-group-message":
				if (msg.status == "success")
				{
					err($("#groupemailError"),"Email message successfully sent.","success");
				}
				else
				{
					err($("#groupemailError"),"Message not sent.","danger");
				}
				break;
				
			/**
			 * A file's public status was updated - we dynamically add an entry or remove its entry 
			 * from users' public file tables
			 */
			case "file-public-change":
				if (msg.status == "success")
                {
					//this was our update
					if (msg.dd.owner == this.getUser())
					{
						$( "#file-manager-dlg" ).modal("hide");
						var dd = this.dd_private_data[msg.dd.id];
						
						//initially public, remove from public listing
						if (dd.ispublic == this.PRIVACY_PUBLIC && msg.dd.ispublic != this.PRIVACY_PUBLIC)
						{
							this.dd_public_data[msg.dd.id].deleteme();
							delete this.dd_public_data[msg.dd.id];
						}
						//show in public data list 
						else if (msg.dd.ispublic == this.PRIVACY_PUBLIC) 
							this.pfiletable.row.add(dd.asArray()).draw( false );

						//update entry in private files
						dd.ispublic = msg.dd.ispublic;
						dd.update();
					}
					//this is a broadcast message
					else
					{
						var id = this._viewingFile(msg.dd.id);
						//check if we are viewing the file - and remove it if necessary
						if (id != -1)
						{
							if (msg.dd.ispublic == this.PRIVACY_ONLYME)
							{
								b.nn_removeTier({"index": id});
								showAlert("File Permissions", "Access to "+msg.dd.trackname+" was revoked.");
							}
							//current user not in a session
							if (this.active_session == undefined)
							{
								if (msg.dd.ispublic == this.PRIVACY_PRIVATE)
								{
									b.nn_removeTier({"index": id});
									showAlert("File Permissions", "Access to "+msg.dd.trackname+" was revoked.");
								}
							}
						}
						
						//check if we are in the same collaborative session
						if (msg.uuid != undefined)
						{
							if (msg.uuid == this.active_session)
							{
								//add the file to private listing
								if (this.dd_private_data[msg.dd.id] == undefined)
								{
									//file is not listed - add it if its NOT ONLYME
									if (msg.dd.ispublic != this.PRIVACY_ONLYME)
									{
										this.filetable.row.add([ msg.dd.id, msg.dd.owner, msg.dd.trackname, msg.dd.description, msg.dd.path,
										                         msg.dd.format, msg.dd.ispublic, msg.dd.indexedby==null?"N/A":msg.dd.indexedby, 
										                         msg.dd.isindex, msg.dd.remote ])
										                         .draw( false );
									}
								}
								else
								{
									//file is listed - delete if its ONLYME
									if (msg.dd.ispublic == this.PRIVACY_ONLYME)
									{
										this.dd_private_data[msg.dd.id].deleteme();
										delete this.dd_private_data[msg.dd.id];
									}
									else
									{
										//update the dd
										this.dd_private_data[msg.dd.id].ispublic = msg.dd.ispublic;
										this.dd_private_data[msg.dd.id].update();
									}
								}
							}
						}
						
						//not already listed
						if (this.dd_public_data[msg.dd.id] == undefined)
						{
							if (msg.dd.ispublic == this.PRIVACY_PUBLIC)
							{
								//add this as a new DD
								this.pfiletable.row.add([ msg.dd.id, msg.dd.owner, msg.dd.trackname, msg.dd.description, msg.dd.path,
								                          msg.dd.format, msg.dd.ispublic, msg.dd.indexedby==null?"N/A":msg.dd.indexedby, 
								                          msg.dd.isindex, msg.dd.remote ])
								                          .draw( false );
							}
						}
						else if (this.dd_public_data[msg.dd.id] != undefined)
						{
							if (msg.dd.ispublic != this.PRIVACY_PUBLIC)
							{
								this.dd_public_data[msg.dd.id].deleteme();
								delete this.dd_public_data[msg.dd.id];
							}
						}
					}
				}
				else
				{
					err($("#fmdlg_response"),msg.error,"warning");
				}
				break;
				
			/**
			 * Broadcast message: a new public access file was uploaded
			 */
			case "new-public-file":
				//add new row to public file table
				this.pfiletable.row.add([ msg.id, msg.owner, msg.trackname, msg.description, msg.path,
            		                      msg.format, msg.ispublic, msg.indexedby==null?"N/A":msg.indexedby, 
            		                      msg.isindex, msg.remote ])
            		                      .draw( false );
				break;

			/**
			 * A change has occurred in the collaborative session and the server is notifying
			 * all the clients to update their statuses
			 */
			case "update-session-users":
				//tier update
				if (msg.tiers !== undefined)
				{
					//don't self update, that's stupid
					if (msg.caller != this.getUser())
					{
						this.synchGenomes(msg.genome, msg.tiers, msg["custom_genome"]);
					}
				}
				//chat relay message
				else if (msg.smsg !== undefined)
				{
					var schat = document.getElementById('schatArea');
					
					var fmsg = "";
					if (msg.raw == true) fmsg = msg.smsg;
					else fmsg = msg.username+"("+msg.nickname+"): "+msg.smsg;
					
					fmsg = "["+msg.time+"] - "+fmsg;
					
					schat.appendChild(document.createTextNode(fmsg));
					schat.appendChild(document.createElement("br"));
					
					var panel = $("#sessionChatPanel"); 
					var indicator = panel.find(".ni");
					if (panel.find(".collapsed").length > 0)
					{
						indicator.text(parseInt(indicator.text())+1);
						indicator.show(100);
					}
					
					var js = $(schat);
				    if(js.length)
				    	js.scrollTop(js[0].scrollHeight - js.height());
				}
				//position/highlight update
				else
				{
					if (msg.username != this.getUser())
					{
						if (msg.highlights.length < b.highlights.length)
						{
							b.clearHighlights();
						}
						
						//recreate set of highlighted regions
						for (var i=0; i<msg.highlights.length; i++)
						{
							b.highlightRegion(msg.highlights[i].chr,msg.highlights[i].min,msg.highlights[i].max);
						}
					}
					var x = [msg];
					this.updateUserList(x);
				}
				break;
				
			/**
			 * A collaborative session options change has occurred
			 */
			case "update-session-options":
				this.processOptions(msg.options);
				break;

			/**
			 * The user has attempted to register a new account
			 */
            case "register":
                if (msg.status == "success")
                {
                    $("#registerDLG").modal('hide');
                    showAlert("Account Activation",
                    		"A message has been sent to the provided email address with activation instructions.");
                }
                else if (msg.status == "fail")
                {
                	err($(".validateTips"),msg.error,"danger");
                }
                break;
                
            /**
             * The user has attmpted to login to an existing account
             */
            case "login":
            	if (msg.status == "success")
                {
            		$('#loginDLG').modal('hide');
            		
            		//set login cookies
            		this.setCookie("username", msg.username);
            		this.setCookie("uauth", msg.uauth);
            		
            		$( "#logoutButton" ).show();
            		$( "#loginButton" ).hide();
            		$( "#registerButton" ).hide();
            		
            		$("#profilemenu").removeClass("hidden");
            		$(".loggedInUser").html(msg.username);
            		
            		//build group list
            		if (msg.groups != undefined)
        			{
            			for (var name in msg.groups)
            			{
            				this.grouplisttable.row.add([name,msg.groups[name].count,
            				                             msg.groups[name].file_count==null?0:msg.groups[name].file_count,
            				                            		 msg.groups[name].access_level,msg.groups[name].access_level,
            				                            		 msg.groups[name].access_level,msg.groups[name].access_level]);
            			}
            			this.grouplisttable.draw();
        			}

            		//get the user's private file listings
            		this.getUserFiles();
            		
            		//process saved configs for the user
            		if (msg["saved_config"] != null)
        			{
            			var conf = JSON.parse(msg["saved_config"]);
            			if (conf != undefined)
        				{
            				for (var name in conf)
        					{
            					this.savedConfigTable.row.add([conf[name].name,conf[name].desc,
            					                               conf[name].genome, conf[name].time, 
            					                               conf[name].tiers]);
        					}
            				this.savedConfigTable.draw();
        				}
        			}
            		
            		//enables/disables buttons requiring login
            		toggleLoginInterfaces(true);
            		
            		//load comments again - user might have left private comments
            		this.loadComments();
                }
            	else if (msg.status == "fail")
                {
                	err($("#loginError"),msg.error,"danger");
                }
            	else if (msg.status == "resumefail")
                {
                	this.deleteCookie("username");
                	this.deleteCookie("uauth");
                	
                	this.me = undefined;
                	
                	$(".loggedInUser").html("");
                	
                	showAlert("Error","Unable to resume your session. Possibly the server has been restarted.");
                }
            	$("#loginform")[0].reset();
                break;

            /**
             * The client has successfully opened a new connection to the server via the websocket
             * Only happens on client initialisation
             */
            case "new-connection":
            	if (this.getUAuth() != undefined && this.getUser() != undefined)
        		{
            		//attempt to restore users session
            		var p = {
        	            "username" : this.getUser(),
        	            "uauth" : this.getUAuth()
        	        };

        			var message = {
        				type: "resume-session",
        	            params: p
        			};
        			this.sendMessage(JSON.stringify(message));
        		}
            	break;
            
            /**
             * The calling user has successfully logged out of their account
             */
            case "logout":
            	this.me = undefined;

            	if (this.active_session)
            		this.leaveSession();
            	
            	this.deleteCookie("username");
            	this.deleteCookie("uauth");
            	
            	$( "#logoutButton" ).hide();
            	$( "#loginButton" ).show();
            	$( "#registerButton" ).show();
                
            	$("#profilemenu").addClass("hidden");
            	$("#sessionPanel").addClass("hidden",250);
            	$(".loggedInUser").html("");
            	
            	this.group_mem_list = {};
            	
            	//erase all groups
            	for (var name in this.grouplist)
        		{
            		this.grouplist[name].deleteme();
            		delete this.grouplist[name];
        		}
            	this.active_group = undefined;
            	
            	//erase all private data
            	for (var dd in this.dd_private_data)
        		{
            		this.dd_private_data[dd].deleteme();
            		delete this.dd_private_data[dd];
        		}
            	
            	//hide login only interfaces
            	toggleLoginInterfaces(false); 
            	$("#sessionChatPanel").addClass("hidden");
            	$("#newComment").attr("disabled","");
            	
            	//reset group view
            	$("#groupoverviewdiv").removeClass("hidden");
    			$("#groupeditordiv").addClass("hidden");
            	break;
                
            /**
             * Retrieved list of files uploaded by the calling user - both public and private
             */
            case "user-files":
            	for (var i=0; i<msg.files.length; i++)
        		{
            		//add new row to the private file table
            		if (this.dd_private_data[msg.files[i].id] == undefined)
            			this.filetable.row.add([ msg.files[i].id, msg.files[i].owner, msg.files[i].trackname, msg.files[i].description, msg.files[i].path,
            		                         msg.files[i].format, msg.files[i].ispublic, msg.files[i].indexedby, msg.files[i].isindex, msg.files[i].remote ]);
        		}
            	this.filetable.draw(false);
            	break;
            	
            /**
             * Retrieved list of public files uploaded by all users of PBrowse
             */
            case "public-files":
            	for (var i=0; i<msg.files.length; i++)
        		{
            		//add new row to the public file table
            		if (this.dd_public_data[msg.files[i].id] == undefined)
            			this.pfiletable.row.add([ msg.files[i].id, msg.files[i].owner, msg.files[i].trackname, msg.files[i].description, msg.files[i].path,
            		                         msg.files[i].format, msg.files[i].ispublic, msg.files[i].indexedby, msg.files[i].isindex, msg.files[i].remote ]);
        		}
            	this.pfiletable.draw();
            	break;
            	
            /**
             * Broadcast message sent to all users viewing the affected file, notifying them that a
             * new comment has been created
             */
            case "make-comment":
            	if (msg.status == "success")
        		{
            		this.updateCommentList([msg.comment],true);
        		}
            	break;
            	
            /**
             * The retrieved list of comments for a particular track
             */
            case "all-track-comments":
            	this.updateCommentList(msg.comments);
            	break;
            	
            /**
             * Broadcast message sent to all users viewing the affected file, notifying them that a
             * comment has been removed by the original author
             */
            case "delete-comment":
            	this._removeComment(this.comments[msg.id]);
            	break;
            	
            /**
             * The calling user has attempted to create a new group
             */
            case "create-group":
            	if (msg.status == "success")
        		{
            		this.grouplisttable.row.add([msg.groupname,1,0,15,15,15,15]).draw(false);
            		err($("#createGroupError"),"Group created.","success");
        		}
            	else
        		{
            		err($("#createGroupError"),"The provided group name is already taken.","warning");
        		}
            	break;
            	
        	/**
             * The user was added to a group by another user
             */
            case "added-to-group":
            	if (msg.status == "success")
        		{
            		//add new row to group list table
            		this.grouplisttable.row.add([msg.groupname,msg.usercount,
            		                             msg.filecount==null?0:msg.filecount,
    		                            		 msg.access_level,msg.access_level,
    		                            		 msg.access_level,msg.access_level])
        			                             .draw(false);
            		
            		showAlert("New Group!","You have been invited to join the group: '"+msg.groupname+"'. " +
            				"View it in the group management page.");
        		}
            	break;
            	
            /**
             * Retrieved list of all the groups the calling user belongs to. The list is stored by
             * the client and updated dynamically
             */
            case "user-groups":
            	for (var name in msg.groups)
        		{
            		if (this.grouplist[name] != undefined)
        			{
            			this.grouplist[name].groupname = name;
            			this.grouplist[name].num_members = msg.groups[name].count;
            			this.grouplist[name].num_files = msg.groups[name].file_count;
            			this.grouplist[name].access_level = msg.groups[name].access_level;
            			
            			this.grouplist[name].update();
        			}
            		else
        			{
            			this.grouplisttable.row.add([name,msg.groups[name].count,
            			                             msg.groups[name].file_count==null?0:msg.groups[name].file_count,
            			                             msg.groups[name].access_level,msg.groups[name].access_level,
            			                             msg.groups[name].access_level,msg.groups[name].access_level])
            			                             .draw(false);
        			}
        		}
            	break;
            	
            /**
             * Update to the permissions of the specified user. The client updates their record of the
             * user's permissions as well
             */
            case "update-group-user-acl":
            	if (msg.status == "success")
        		{
            		if (this.group_mem_list[msg.groupname] != undefined)
            		{
            			if (this.group_mem_list[msg.groupname][msg.username] != undefined)
            			{
            				var mem = this.group_mem_list[msg.groupname][msg.username];
            				mem.setacl(msg.acl);
            			}
            		}
        		}
            	
            	if (msg.username == this.getUser())
        		{
            		this.grouplist[msg["groupname"]].access_level = msg.acl;
            		this.grouplist[msg["groupname"]].update();
            		
            		if (this.active_group == msg.groupname)
        			{
            			this.grouplist[msg["groupname"]].hideMenus();
            			$('.nav-tabs a[href="#groupmembers"]').tab('show');
        			}
        		}
            	break;
            	
            /**
             * The user has attempted to share a file with the specified group. The client maintains
             * and dynamically updates the list with received updates.
             */
            case "share-group-file":
            	if (msg.status == "success")
        		{
            		this.grouplist[msg["groupname"]].num_files += 1;
            		this.grouplist[msg["groupname"]].update();
            		
            		//invoked by other user -- notify current user
            		if (msg.caller != this.getUser())
            		{
            			showAlert("Shared Files","User: "+msg.caller+" has shared a new file with '"+msg.groupname+"', with track-name: "+msg.file.trackname+". Add the file " +
            			"as a track by visiting the <a href='#' data-toggle='modal' data-target='#groupmanagerDLG'>group management page.</a>");
            		}
            		else
        			{
            			$( "#file-manager-dlg" ).modal("hide");
        			}
            		
            		//the user has not loaded the group files beforehand
            		if (this.dd_group_data[msg.groupname] == undefined)
            			break;
            		
            		this.gfiletable.row.add([ msg.file.id, msg.file.owner, msg.file.trackname, 
        			                          msg.file.description, msg.file.path,
	            		                      msg.file.format, msg.file.ispublic, 
	            		                      msg.file.indexedby, msg.file.isindex, msg.file.remote, msg.groupname ])
	            		                      .draw( false );
            		
            		//invoke hidden group filter
            		if (this.active_group != undefined)
            			this.gfiletable.column(10).search(this.active_group);
        		}
            	else
        		{
            		err($("#fmdlg_response"),"The file could not be shared with the selected group.","info");
        		}
            	break;
            	
            /**
             * A new user was added to the group, the client maintains a list of users and updates it
             * dynamically
             */
            case "add-group-user":
            	if (msg.status == "success")
        		{
            		if (this.group_mem_list[msg["groupname"]] != undefined)
        			{
            			this.groupmembertable.row.add([msg["username"],msg["groupname"],msg["access_level"],
            			                               msg["access_level"],msg["access_level"],msg["access_level"],
            			                               msg["access_level"]])
            			                               .draw(false);
            			
            			this.grouplist[msg["groupname"]].num_members += 1;
            			this.grouplist[msg["groupname"]].update();
            			
            			err($("#addNewUserError"),"The user has been successfully added.","success");
        			}
        		}
            	else
        		{
            		err($("#addNewUserError"),"The user does not exist. Or already belongs to the group.","warning");
        		}
            	break;
            	
            /**
             * Retrieved list of all the files shared with a group and current group members.
             * The client dynamically updates this list when changes occur
             */
            case "group-info":
            	if (this.dd_group_data[msg.groupname] == undefined)
            	{
            		this.dd_group_data[msg.groupname] = {};
            	}
            	
            	for (var i=0; i<msg.files.length; i++)
        		{
            		if (this.dd_group_data[msg.groupname][msg.files[i].id] == undefined)
        			{
            			this.gfiletable.row.add([ msg.files[i].id, msg.files[i].owner, msg.files[i].trackname, 
            			                          msg.files[i].description, msg.files[i].path,
		            		                      msg.files[i].format, msg.files[i].ispublic, 
		            		                      msg.files[i].indexedby, msg.files[i].isindex, msg.files[i].remote, msg.groupname ])
		            		                      .draw( false );
        			}
    			}
            	
            	for (var i=0; i<msg.users.length; i++)
        		{
            		var mem = null;
            		if (this.group_mem_list[msg.groupname] != undefined)
        			{
            			if (this.group_mem_list[msg.groupname][msg.users[i]["username"]] != undefined)
            			{
                			mem = this.group_mem_list[msg.groupname][msg.users[i]["username"]];
                			mem.setacl(msg.users[i]["access_level"]);
            			}
        			}
            		else this.group_mem_list[msg.groupname] = {};
            			
            		if (mem == null)
        			{ 
            			this.groupmembertable.row.add([msg.users[i]["username"],msg.groupname,msg.users[i]["access_level"],
            			                               msg.users[i]["access_level"],msg.users[i]["access_level"],
            			                               msg.users[i]["access_level"],msg.users[i]["access_level"]])
            			                               .draw(false);
        			}
            	}
            	break;
            	
            /**
             * The group was deleted by the owner, all the members are removed and all shared files are no
             * longer accessible - though they are still accessible to the original uploaders
             */
            case "group-deleted":
            	if (msg.status == "success")
        		{
            		this.grouplist[msg.groupname].deleteme();
            		delete this.grouplist[msg.groupname];
            		
            		if (this.active_group == msg.groupname)
        			{
            			this.active_group = undefined;
            			
            			if ($("#groupManagerPanel").css("display") != "none")
        				{
            				//return to overview screen
            				$("#groupoverviewdiv").removeClass("hidden");
            				$("#groupeditordiv").addClass("hidden");
        				}
        			}
            		
            		if (msg.caller != this.getUser())
        			{
            			showAlert("Group Management","The '"+msg.groupname+"' group was deleted by the owner.");
        			}
        		}
            	else
        		{
            		//should only fail on illegal invocation, i.e. by illegitimate user
        		}
            	break;
            	
            /**
             * The specified user was removed from the group, the client dynamically updates this list
             */
            case "remove-group-user":
            	if (this.group_mem_list[msg.groupname] != undefined)
    			{
        			if (this.group_mem_list[msg.groupname][msg.username] != undefined)
        			{
            			this.group_mem_list[msg.groupname][msg.username].deleteme();
            			delete this.group_mem_list[msg.groupname][msg.username]; 
            			
            			this.grouplist[msg.groupname].num_members -= 1;
            			this.grouplist[msg.groupname].update();
        			}
    			}
            	
            	//user leaving group case
            	if (this.getUser() == msg.username)
        		{
            		this.grouplist[msg.groupname].deleteme();
            		delete this.grouplist[msg.groupname];
            		
            		$("#groupoverviewdiv").removeClass("hidden");
        			$("#groupeditordiv").addClass("hidden");
        		}
            	break;
            	
            /**
             * The access requirement for a particular file belonging to a specified group was modifed
             * NOW DEPRECATED
             */
            case "update-group-file-acl":
            	var dd = this.group_shared_files[msg.groupname][msg.fileid];
            	
            	//something went wrong - the file must exist before being updated
            	if (dd == undefined) return;
            	
            	//update the entry
            	dd.setACL(msg.acl);
            	break;
            	
            /**
             * The file shared with a particular group was un-shared
             */
            case "remove-group-file-acl":
            	this.grouplist[msg.groupname].num_files -= 1;
            	this.grouplist[msg.groupname].update();

            	//the user has not loaded the group files beforehand
            	if (this.dd_group_data[msg.groupname] == undefined)
            		break;
            	
            	var dd = this.dd_group_data[msg.groupname][msg.fileid];
            	
            	//something went wrong - the file must exist before being updated
            	if (dd == undefined) break;
            	
            	//update the entry
            	dd.deleteme();
            	delete this.dd_group_data[msg.groupname][msg.fileid];
            	
            	break;
        	
            /**
             * The user has changed their password
             */
            case "change-password":
            	if (msg.status == "success")
        		{
            		err($("#changePasswordError"),"Password SUCCESSFULLY changed!","success");
        		}
            	else
        		{
            		err($("#changePasswordError"),"Password was NOT changed!","danger");
        		}
            	break;
            	
        	/**
        	 * The user has invoked a request to change their recovery email address to a new one
        	 */
            case "change-email":
            	if (msg.status == "success")
        		{
            		err($("#changeEmailError"),"An email was sent to the entered address with instructions.","success");
        		}
            	else
        		{
            		err($("#changeEmailError"),"The operation failed!","danger");
        		}
            	break;
            	
        	/**
        	 * The user changed their nickname to the specified one
        	 */
            case "change-nickname":
            	if (msg.status == "success")
        		{
            		err($("#changeNicknameError"),"Nickname was changed!","success");
            		if (this.session_users[msg.username] != undefined)
        			{
            			this.session_users[msg.username].nickname = msg.nickname;
            			this.session_users[msg.username].update();
            			
            			this._updateStatus();
        			}
        		}
            	else
        		{
            		err($("#changeNicknameError"),"Password was NOT changed!","danger");
        		}
            	break;
            	
        	/**
        	 * The user requested a password reset and an email was sent to the corresponding address
        	 */
            case "reset-password":
            	if (msg.status == "success")
        		{
            		$('#loginDLG').modal('hide');
            		showAlert("Password Reset Request","An email has been sent to the associated address " +
            				"with reset instructions.");
        		}
            	else
        		{
            		showAlert("Password Reset Request","Invalid username or email address given.");
        		}
            	break;
            	
        	/**
        	 * this message only confirms the success of the nomination, an user update
        	 * immediately follows it which formalises the change for all session users
        	 */
            case "session-nominate-leader":
            	if (msg.status == "success")
        		{
            		if (msg.caller == this.getUser())
        			{
            			//notification text
            			this.sendSessionMsg(msg.caller+" has nominated "+msg.leader+" as the new leader.", true);
            			$("#sum-dlg").modal('hide');
        			}
            		
            		this.updateUserList(msg["user-status"]);
        		}
            	break;
            	
        	/**
        	 * The session leader blacklisted a particular user, if that user was in the current session,
        	 * they are kicked from it and prevented from returning. If the leader removes them from the 
        	 * blacklist, they may rejoin
        	 */
            case "session-blacklist-update":
            	if (msg.status == "success")
        		{
            		if (this._amILeader())
        			{
            			$("#sum-dlg").modal('hide');
            			if (msg.option == "deny")
        				{
            				err($("#sessionBlacklistUserError"),"The user has been blacklisted!","success");
            				//notification text
            				this.sendSessionMsg(this.getUser()+" has blacklisted "+msg.user, true);
        				}
            			else if (msg.option == "allow")
        				{
            				this.sendSessionMsg(this.getUser()+" has whitelisted "+msg.user, true);            				
            				err($("#sessionBlacklistUserError"),"The user has been removed from the blacklist!","success");
        				}
        			}
            		if (msg.option == "allow")
        			{
            			if (this.session_blacklist[msg.user] != undefined)
            				delete this.session_blacklist[msg.user];
        			}
            		else if (msg.option == "deny")
        			{
            			this.session_blacklist[msg.user] = 0;
        			}
            		this.drawBlacklist();
        		}
            	else
        		{
            		err($("#sessionBlacklistUserError"),"The user could not be blacklisted!","warning");
        		}
            	break;
            	
        	/**
        	 * The session leader invited a specific user to join their collaborative session
        	 */
            case "session-invite-user":
            	if (msg.status == "success")
        		{
            		if (msg.user == this.getUser())
        			{
            			$("#joinSessionDLG").modal("show");
            			$("#sessionIDText").val(msg.uuid);
            			$("#joinSessionCode").val(msg.code);
            			
            			err($("#joinsessionError"),"You have been invited to join this sesion. Press 'Join Session' to continue.","info");
        			}
            		else
        			{
            			err($("#sessionInviteUserError"),"Invitation was sent!","success");
        			}
        		}
            	else
        		{
            		err($("#sessionInviteUserError"),msg.error,"warning");
        		}
            	break;
            	
            /**
             * A particular user was kicked by the leader of a session
             */
            case "session-kick-user":
            	if (msg.status == "success")
        		{
            		if (msg.user == this.getUser())
        			{
            			if (msg.option == "deny")
        				{
            				showAlert("Blacklisted!","You have been blacklisted from the session!");
        				}
            			else
        				{
            				showAlert("Kicked!","You have been kicked from the session!");
            				this.sendSessionMsg(this.getUser()+" has been kicked.", true);
        				}
            			this.leaveSession();
        			}
            		else
        			{
            			$("#sum-dlg").modal('hide');
        			}
        		}
            	break;
            	
            case "update-saved-config":
            	if (msg.status == "success")
        		{
            		err($("#cfg-save-error"),"Configuration saved successfully.","success");
        		}
            	else
        		{
            		err($("#cfg-save-error"),"Configuration could not be saved!","danger");
        		}
            	break;
            	
        	/**
        	 * Unimplemented message - should never occur
        	 */
            default:
				showAlert("Error","DEFAULT: "+event.data);
				break;
		}
	};
	/**
	 * Sets the new cookie given by the name and value parameters, to expire automatically
	 * after 24 hours
	 */
	this.setCookie = function (name, value)
	{
		var d = new Date();
	    d.setTime(d.getTime() + (24*60*60*1000));
	    var expires = "expires="+d.toUTCString();
		document.cookie = name+"="+value+";"+expires;
	};
	/**
	 * Deletes the cookie indentified by the name parameter, setting its expiry date to one
	 * in the past
	 */
	this.deleteCookie = function(name)
	{
		document.cookie = name+"=; expires=Thu, 01 Jan 1970 00:00:00 UTC";
	};
	/**
	 * Updates the comment tree with the values provided in the 'data' parameter. if the 'flag'
	 * parameter is set, indicate that the comments are new and notify the user 
	 */
	this.updateCommentList = function (data, flag)
	{
		var comments = this.comments;

		//rootnode -- global for testing
		var root = this.rootComment;
		if (root == null)
		{
			root = new CommentNode(0,0,0);
		}
				
		for (var c in data)
		{
			var com = data[c];
			
			//selectively ignore comment updates if we don't view the relevant file
			if (this.tierIDMap[com.trackref] == undefined)
			{
				continue;
			}
			
			//attach trackref node to root
			if (root.children[com.trackref] == undefined)
			{
				var trefnode = new CommentNode(0,com.trackref,0);
				root.addc(com.trackref,trefnode);
				trefnode.parent = root;
			}
			
			//add chr node to the tracknode
			var chr_node = root.children[com.trackref].children[com.chr];
			if (chr_node == undefined)
			{
				chr_node = new CommentNode(1,0,com.chr);
				root.children[com.trackref].addc(com.chr,chr_node);
				chr_node.parent = root.children[com.trackref];
				
				if (chr_node.parent.expanded)
				{
					//make visible
					chr_node.me.show();
				}
			}
			
			if (comments[com.id] == undefined)
			{
				//this comment is new
				comments[com.id] = new Comment(com.id, com.trackref, com.chr, com.startpos, com.endpos, com.creator, com.ctext, com.ispublic);
				//add the comment into the tree under its chr node mapped by comment id
				chr_node.addc(com.id,comments[com.id]);
				comments[com.id].parent = chr_node;
				
				//new comments are only generated by 'make comment'
				if (flag != undefined)
					comments[com.id].setnew(1);
			}
			
			if (comments[com.id].parent.expanded)
			{
				comments[com.id].me.show();
			}
		}
		
		//populate comment list field
		if (!root.indom)
		{
			root.indom = true;
			$("#commentlist").append(root.me);
		}
		this.rootComment = root;
	};
	/**
	 * Function handling all changes to the user list within a collaborative session - takes an
	 * array of status data as input and updates the user list with new values
	 */
	this.updateUserList = function (data)
	{
		if (data == undefined)
		{
			return;
		}
		
		for (var i=0; i<data.length; i++)
		{
			var user = this.session_users[data[i].username];
			
			var isl = data[i]["isLeader"]!==undefined?data[i]["isLeader"]:true;
			
			//if new user - initiate new IUser
			if (user === undefined)
			{
				user = new IUser(data[i]["username"], 
						data[i]["nickname"], 
						isl,
						data[i]["chr"],
						data[i]["start_win"], 
						data[i]["end_win"]);
				
			}
			//otherwise modify him
			else
			{
				user.chr = data[i]["chr"];
				user.start = data[i]["start_win"];
				user.end = data[i]["end_win"];
				user.isleader = isl;
			}
			
			//if this user has become leader, enable admin panel
			if (isl == true && user.username == this.getUser())
			{
				//enable admin panel
				$("#adminControllerDiv").show();
				$("#destroySessionBtn").removeAttr("disabled");
				
				$("#sessionLeaderIndicator").text("| Leader");
			}
			else if (!this._amILeader())
			{
				$("#adminControllerDiv").hide();
				$("#destroySessionBtn").attr("disabled","");
				
				$("#sessionLeaderIndicator").text("");
			}
			
			//invoke update on user html
			user.update();
			
			//if user isn't in the user-list DOM already, insert him
			if (!user.indom)
			{
				$("#userListDiv").append(user.me);
				user.indom = true;
			}

			//store new user into user-list
			this.session_users[data[i].username] = user;

			if (this._follower_sync == false)
			{
				//dont sync if option set
				continue;
			}

			//subscribed user updates followers
			if ( user.subscribed == true )
			{
				//only sync update if it was invoked by a human
				if ( user.username == this.getUser() || !data[i].human )
					continue;
				
				if ((b.viewEnd|0) == user.end && (b.viewStart|0) == user.start)
					continue;
				
				//no notify
				b.nn_setLocation(user.chr,user.start,user.end);
				this._updateStatus();
			}
		}
	};
	/**
	 * The server encountered a problem and closed the websocket connection for this user
	 * or the client navigated away from the page
	 */
	this.onclose = function (event)
	{
		console.warn("SOCKET CLOSED:");
		console.warn(JSON.stringify(event));
		if (event.reason.length > 0)
		{
			showAlert("Error","The server has closed the websocket connection. Reconnecting automatically in 5 seconds.");
		}
		window.setTimeout(function(){
			location.reload();
		}, 5000);
	};
	
	/**
	 * Changes the viewed genome to the one specified by 'key' parameter
	 */
	this.changeGenome = function(key, callback)
	{
		this.active_genome = key;
		
		var config = genomes[key];
		config["maxWorkers"] = 2;
		
		b=new Browser(config);
		b.addInitListener(this.onbrowserinit.bind(this));
		
		if (callback != undefined)
		{
			b.addInitListener(callback);
		}
	};
	
	//dalliance browser init listener, initialises several listener functions
	this.onbrowserinit = function()
	{
		var _t = this;
		
		b.addTierListener(_t.onTierChange.bind(this));
		b.addViewListener(_t.onViewChange.bind(this));
		b.addTierSelectionListener(function (sel) {
			if (sel.length == 1)
			{
				if (b.tiers[sel[0]].dasSource._id != undefined)
				{
					if (_t.getUser() != undefined)
						$(".tierCommentButton").parent().removeAttr("disabled");
				}
				else
				{
					$(".tierCommentButton").parent().attr("disabled","");
				}
			}
		});
		b.storeTierStatus();
		_t._updateStatus();
		
		_t.loadComments();
		this._numloadedtracks = b.tiers.length;
	};
	
	//tier change listener, called whenever the user initiates a browser track modification
	this.onTierChange = function(e)
	{
		var _t = this;
		
		if (this._numloadedtracks != b.tiers.length)
		{
			_t.loadComments();
			this._numloadedtracks = b.tiers.length;
		}
		
		 //only trigger update when this is the last request received in the past 200ms
		if (_t.tiers_timeout == null)
		{
			_t.tiers_timeout = setTimeout(_t.updateTiers, 200);
		}
		else
		{
			window.clearTimeout(_t.tiers_timeout);
			_t.tiers_timeout = setTimeout(_t.updateTiers, 200);
		}
	};
	
	//view change listener, called when the user moves the view via mouse, keyboard, or other interface
	this.onViewChange = function(chr, min, max)
	{
		var _t = this;
		
	    //only trigger update when this is the last request received in the past 750ms
		if (_t.update_timeout == null)
		{
			_t.update_timeout = setTimeout(function(){_t.updateStatus();}, 750);
		}
		else
		{
			clearTimeout(_t.update_timeout);
			_t.update_timeout = setTimeout(function(){_t.updateStatus();}, 750);
		}
	};
	
	this.init = function()
	{
		//init the browser
		cm.changeGenome(cm.active_genome);
		
		//ROW processing functions for all datatables used in pbrowse
		var rowcreate = function ( row, data, index ) 
        {
        	var desc = $('td', row).eq(3);
        	var indexid = $('td', row).eq(7);
        	$('td', row).eq(2).css("word-break", "break-all"); //trackname
        	$('td', row).eq(4).css("word-break", "break-all"); //path
        	desc.css("word-break", "break-all");

        	//remotely located file - process metadata in format
        	if (data[9] != undefined)
        	{
                if (data[9] == true)
                {
                    var meta = JSON.parse(data[5]);
                    $('td', row).eq(5).html(meta.data_format);

                    //modify path into normal /studyid/file format
                    $('td', row).eq(4).html("/"+meta.studyid+"/"+meta.filename);
                }
    		}
        	
        	var ispublic = $('td', row).eq(6);
        	if (data[6] == 4) // PUBLIC
        		ispublic.html(
        		'<span class="glyphicon glyphicon-globe p_public" title="Everyone" aria-hidden="true"></span>');
        	else if (data[6] == 2) // PRIVATE
        		ispublic.html(
        		'<span class="glyphicon glyphicon-eye-open p_private" title="Yourself and Collaborators" aria-hidden="true"></span>');
        	else if (data[6] == 1) // ONLYME
        		ispublic.html(
        		'<span class="glyphicon glyphicon-eye-close p_onlyme" title="Only Yourself" aria-hidden="true"></span>');
        	
        	if (desc.html().indexOf("#") != -1)
        	{
            	var t = desc.html().substr(1).split(":");
            	var tag = "<span class='label label-"+t[0].toLowerCase()+"'>"+t[0]+"</span>";
            	desc.html(tag+(t[1]==undefined?"":t[1]));
        	}
        	if (indexid.html()=="0")
        		indexid.html("N/A");
        };
        var privatedata = function( row, data, index )
        {
        	rowcreate( row, data, index );
        	var dd = new DataDesc(data, row, this.filetable, 0);
        	this.dd_private_data[data[0]] = dd;
        }.bind(this);
        var publicdata = function( row, data, index )
        {
        	rowcreate( row, data, index );
        	var dd = new DataDesc(data, row, this.pfiletable, 1);
        	this.dd_public_data[data[0]] = dd;
        }.bind(this);
        var groupdata = function( row, data, index )
        {
        	rowcreate( row, data, index );
        	var dd = new DataDesc(data, row, this.gfiletable, 2);
        	
        	if (this.dd_group_data[data[10]] == undefined)
        		this.dd_group_data[data[10]] = {};
        	
        	this.dd_group_data[data[10]][data[0]] = dd;
        }.bind(this);
		
        //list of all user's uploaded files
		this.filetable = $("#fileTable").DataTable( {
			"autoWidth": false,
	        "createdRow": privatedata,
	        select: true,
	        "columnDefs": [
	            { className: "dt_hidden_column", "targets": [7,8,9] }
	        ]
	    });
		//list of the currently publicly available files
		this.pfiletable = $("#publicfileTable").DataTable( {
			"autoWidth": false,
	        "createdRow": publicdata,
	        select: true,
	        "columnDefs": [
	            { className: "dt_hidden_column", "targets": [6,7,8,9] }
	        ]
	    });
		//list of all the files belonging to all the groups the user belongs to
		this.gfiletable = $("#groupFileTable").DataTable( {
			"autoWidth": false,
	        "createdRow": groupdata,
	        select: true,
	        "dom": '<<"row"<"col-sm-6"l><"col-sm-6"f>><t>p>',
	        "columnDefs": [
	            { className: "dt_hidden_column", "targets": [6,7,8,9,10] }
	        ]
	    });
		//list of all groups the user belongs too
		this.grouplisttable = $("#groupRoleTable").DataTable( {
			"autoWidth": false,
	        "createdRow": function ( row, data, index ) 
	        {
	        	//add row classes
	        	$('td', row).eq(0).addClass("grow_gn");
	        	$('td', row).eq(1).addClass("grow_nm");
	        	$('td', row).eq(2).addClass("grow_nf");
	        	//checkboxes
	        	$('td', row).eq(3).html("<div class='grpn checkbox disabled'><label><input class='grow_p_r' type='checkbox' value='' disabled></label></div>");
	        	$('td', row).eq(4).html("<div class='grpn checkbox disabled'><label><input class='grow_p_mf' type='checkbox' value='' disabled></label></div>");
	        	$('td', row).eq(5).html("<div class='grpn checkbox disabled'><label><input class='grow_p_mu' type='checkbox' value='' disabled></label></div>");
	        	$('td', row).eq(6).html("<div class='grpn checkbox disabled'><label><input class='grow_p_own' type='checkbox' value='' disabled></label></div>");
	        	
    			this.grouplist[data[0]] = new GroupRow(data[0],data[1],data[2],data[3],row,this.grouplisttable);
	        }.bind(this),
	    });
		//handles group member display for all groups
		this.groupmembertable = $("#groupMemberTable").DataTable( {
			"autoWidth": false,
	        "createdRow": function ( row, data, index ) 
	        {
	        	//add row classes
	        	$('td', row).eq(2).addClass("gm_acl");
	        	//checkboxes
	        	$('td', row).eq(3).html("<div class='grpn checkbox disabled'><label><input class='gm_p_r' type='checkbox' value='' disabled></label></div>");
	        	$('td', row).eq(4).html("<div class='grpn checkbox disabled'><label><input class='gm_p_mf' type='checkbox' value='' disabled></label></div>");
	        	$('td', row).eq(5).html("<div class='grpn checkbox disabled'><label><input class='gm_p_mu' type='checkbox' value='' disabled></label></div>");
	        	$('td', row).eq(6).html("<div class='grpn checkbox disabled'><label><input class='gm_p_own' type='checkbox' value='' disabled></label></div>");
	        	
	        	if (this.group_mem_list[data[1]] == undefined) this.group_mem_list[data[1]] = {};
	        	
	        	this.group_mem_list[data[1]][data[0]] = new GroupMember(data[1], data[0], data[2], row, this.groupmembertable);
	        	
	        }.bind(this),
	        "dom": '<<"row"<"col-sm-6"l><"col-sm-6"f>><t><"row"<"col-sm-8 help-block"><"col-sm-4"p>>>',
	        "columnDefs": [
	            { className: "dt_hidden_column", "targets": [1,2] }
	        ]
	    });
		$("#groupmembers").find(".help-block").text("*R = Read Files, *F = Manage Files, *U = Manage Users, *O = Group Owner");
		//handles public sessions display and search
		this.publicsessiontable = $("#publicSessionTable").DataTable( {
			"autoWidth": false,
	        "createdRow": function ( row, data, index ) 
	        {
	        	this.public_sessions_list[data[0]] = new PublicSession(
	        			data[0], data[1], data[2], data[3], row, this.publicsessiontable);
	        }.bind(this),
	    });
		
		this.savedConfigTable = $("#savedConfigTable").DataTable( {
			"autoWidth": false,
	        "createdRow": function ( row, data, index ) 
	        {
	        	//override existing named entries
	        	var ex = this.saved_config[data[0]];
	        	if (ex != undefined)
        		{
	        		ex.deleteme();
	        		delete this.saved_config[data[0]];
        		}
	        	this.saved_config[data[0]] = new SavedConfig(data, row, this.savedConfigTable);
	        }.bind(this),
	        "columnDefs": [
   	            { className: "dt_hidden_column", "targets": [4] }
   	        ]
	    });
		
		this.getPublicFiles();
		this.updateCommentList();
	};
	
	//open the socket, triggering the initialization chain of events
	this.sock = new WebSocket("wss://"+window.location.hostname+":8443/connect");
	//set handlers for socket operation
	this.sock.onopen = this.init.bind(this);
	this.sock.onmessage = this.onmessage.bind(this);
	this.sock.onclose = this.onclose.bind(this);
	this.sock.onerror = function (e){
		console.warn(JSON.stringify(e));
	};
};

//init the collabmanager
var cm = new CollabManager();

//test synchronizer
function Tester()
{
	//create second socket
	this.sock = new WebSocket("wss://"+window.location.hostname+":8443/connect");
	
	this.group = undefined;
	this.question = 0;
	this.isleader = false;
	
	this.setGroup = function (groupname,isleader)
	{
		var p = {
            "groupname" : groupname,
            "isleader" : isleader,
        };
		this.group = groupname;
		this.isleader = isleader;

		var message = {
			type: "TEST-set-group",
            params: p,
		};
		this.sendMessage(JSON.stringify(message));
	};
	
	this.beginQuestion = function (qid)
	{
		var p = {
            "qid" : qid,
        };
		this.group = groupname;

		var message = {
			type: "TEST-begin-question",
            params: p,
		};
		this.sendMessage(JSON.stringify(message));
	};
	
	this.onmessage = function (event)
	{
		var msg = JSON.parse(event.data);
		switch (msg.rtype)
		{
		case "TEST-begin-question":
			currentQuestion = msg.qid;
			$(".startResponse").click();
			break;
		}
	};
	
	this.onclose = function (event)
	{
		showAlert("Tests","Testing completed.");
	};
	
	//init socket connection and callbacks
	this.sock.onmessage = this.onmessage.bind(this);
	this.sock.onclose = this.onclose.bind(this);
};