/**
 * The user instance, stores view related parameters and allows the user-list
 * to be easily updated.
 */
function IUser(username, nickname, isleader, _chr, _start, _end) 
{
    this.chr = _chr;
    this.start = _start;
    this.end = _end;
    this.username = username;
    this.nickname = nickname;
    this.isleader = isleader;
    this.subscribed = false;
   
    this.me = $(
    		"<tr>" +
    		"<td iusr_username></td>" +
    		"<td iusr_pos></td>" +
    		"<td ><span style='padding-top: 1px;' iusr_isleader class='glyphicon glyphicon-star'></span></td>" +
    		"<td><div class='grpn checkbox'><label><input iusr_subscribed type='checkbox' value=''></label></div></td>" +
    		"</tr>"
    		);
    
	this.me.on("click",function(event) 
	{
		var _t = this;
		event.stopPropagation();
		if (cm.getUser() == _t.username)
			return;
		
		//set title
		$("#sum_title").text(_t.username+"("+_t.nickname+")");
		
		//set a new leader
		var sum_sl = $("#sum_setleader");
		sum_sl.off("click");
		sum_sl.on("click",function(e){
			cm.sessionNominateLeader(_t.username);
		});
		//kick the user
		var sum_ku = $("#sum_kickuser");
		sum_ku.off("click");
		sum_ku.on("click",function(e){
			cm.sessionKickUser(_t.username);
		});
		//blacklist the user
		var sum_blu = $("#sum_blacklistuser");
		sum_blu.off("click");
		sum_blu.on("click",function(e){
			cm.sessionBlacklistUpdate(_t.username, "deny");
		});
		//goto the position and synchronise with the user
		var sum_fu = $("#sum_follow_user");
		sum_fu.off("click");
		sum_fu.on("click",function(e){
			cm.session_users[_t.username].subscribed = true;
			_t.e_subscribed[0].checked = true;
			b.setLocation(_t.chr,_t.start,_t.end);
			$("#sum-dlg").modal('hide');
		});
		
		if (!cm._amILeader())
		{
			$("#sum_leaderopts").addClass("hidden");
			$("#sum_followeropts").removeClass("hidden");
		}
		else
		{
			$("#sum_leaderopts").removeClass("hidden");
			$("#sum_followeropts").addClass("hidden");
		}
		
		//simply goto the position
		$("#sum_gotopos").off("click");
		$("#sum_gotopos").on("click",function(e){
			b.setLocation(_t.chr,_t.start,_t.end);
			$("#sum-dlg").modal('hide');
		});
		
		$("#sum-dlg").modal('show');
	}.bind(this));
	
    //if the user is inserted into the user-list already
    this.indom = false;
    
    this.e_username = this.me.find('[iusr_username]');
    this.e_pos = this.me.find('[iusr_pos]');
    this.e_isleader = this.me.find('[iusr_isleader]');
    this.e_subscribed = this.me.find('[iusr_subscribed]');
    
    //if leader default to subscribed on
    if (this.isleader)
	{
    	this.subscribed = true;
    	this.e_subscribed[0].checked = true;
	}
    
    this.e_subscribed.on("click", function(event) 
	{
		event.stopPropagation();
		this.subscribed = this.e_subscribed[0].checked;
	}
    .bind(this));
    
    //update DOM of element
    this.update = function() 
    {
    	this.e_username.text(this.username+"("+this.nickname+")");
    	this.e_pos.text(this.chr+":"+this.start+".."+this.end);
    	
    	if (this.isleader)
		{
    		this.e_isleader.addClass("glyphicon-star");
		}
    	else
    	{
    		this.e_isleader.removeClass("glyphicon-star");
		}
    };
    
    this.deleteme = function() {
    	this.me.remove();
    };
}

/**
 * The object representing rows in the public session interface. Provides means of interacting with
 * the elements and joining the session.
 */
function PublicSession(sid, name, num_users, entrycode, row, parenttable)
{
	this.sid = sid;
	this.name = name;
	this.num_users = num_users;
	this.entrycode = entrycode;
	
	this.me = $(row);
	this._row = row;
	this.parent = parenttable;
	
	this.onclick = function(event)
	{
		event.stopPropagation();
		
		var dlg = $("#public-session-dlg");
		if (this.name.length > 0)
		{
			dlg.find('.modal-title').text(this.name);
		}
		else dlg.find('.modal-title').text("Public Session");
		
		dlg.find("#ps_sid").val(this.sid);
		if (this.entrycode == false)
		{
			dlg.find("#ps_passcode_row").addClass("hidden");
		}
		else dlg.find("#ps_passcode_row").removeClass("hidden");
		
		$("#ps_joinbtn").off("click");
		$("#ps_joinbtn").on("click", function(event) {
			cm.joinSession($("#ps_sid").val(),$("#ps_passcode").val());
			dlg.modal('hide');
		});
		
		dlg.modal('show');
	};
	this.me.on("click",this.onclick.bind(this));
	
	this.deleteme = function()
	{
		this.parent.row(this._row).remove().draw();
	};
}

/**
 * The object representing rows in the public session interface. Provides means of interacting with
 * the elements and joining the session.
 */
function SavedConfig(data, row, parenttable)
{
	this.name 	= data[0];
	this.desc 	= data[1];
	this.genome = data[2];
	this.time 	= data[3];
	this.tiers 	= data[4];
	
	this.me = $(row);
	this._row = row;
	this.parent = parenttable;
	
	this.onclick = function(event)
	{
		var _t = this;
		event.stopPropagation();
		
		var dlg = $("#configManagerDLG");
		
		var loader = dlg.find("#cfgman-load");
		loader.off("click");
		loader.on("click",function(e){
			var _settiers_ = function() {
				b.nn_removeAllTiers();
				
				var tiers = JSON.parse(_t.tiers);
				cm.setTiers(tiers);
				dlg.modal('hide');
				
				b.notifyTier();
			};
			
			if (cm.active_genome != _t.genome)
			{
				if (genomes[_t.genome] == undefined)
				{
					//this is a custom genome
					var gen = genomes["generic"];
					gen.coordSystem.speciesName = _t.genome;
					
					//add each tier as a config object so dalliance loads them on init
					var tiers = JSON.parse(_t.tiers);
					for (var i=0; i<tiers.length; i++)
					{
						gen.sources.push(tiers[i].source);
					}
					
					cm.changeGenome("generic", function(){
						//load our modified generic genome again 
						//so that any style configuration is preserved
						_settiers_();
					});
				}
				else
				{
					cm.changeGenome(_t.genome, function(){
						//after browser init
						_settiers_();
					});
				}
			}
			else
			{
				_settiers_();
			}
		});
		
		var delbtn = dlg.find("#cfgman-delete");
		delbtn.off("click");
		delbtn.on("click",function(e){
			delete cm.saved_config[_t.name];
			_t.deleteme();
			dlg.modal('hide');
			
			//push update
			cm.updateSavedConfig();
		});
		
		dlg.modal('show');
		
		//find existing modal-backdrop
		var modalbd = $(".modal-backdrop:visible");
		modalbd[modalbd.length-1].style.zIndex = dlg.css("z-index")-1;
		
		//keep modal scrolling active
		dlg.on('hidden.bs.modal', function(event) {
			if (modalbd.length >= 1)
				$("body").addClass("modal-open");
		});
	};
	this.me.on("click",this.onclick.bind(this));
	
	this.deleteme = function()
	{
		this.parent.row(this._row).remove().draw();
	};
}

/**
 * The lowest level comment node, contains the actual comment data, position, etc.
 * All comments are grouped into a tree-style hierarchy under CommentNode objects. 
 */
function Comment(id, trackref, chr, start, end, creator, ctext, ispublic) 
{
    this.id = id;
    this.trackref = trackref;
    this.chr = chr;
    this.start = start;
    this.end = end;
    this.creator = creator;
    this.ctext = ctext;
    this.ispublic = ispublic;

    this.parent = null;
    this.setnew = function(t)
    {
    	if (t > 0)
		{
    		this.isnew = true;
    		this.me.addClass("ne");
		}
		else this.isnew = false;
    	
    	if (this.parent != undefined)
		{
    		this.parent.newchildren(t);
		}
    };
    this.isnew = false;
    
    this.me = $("<li class='nc' style='display:none'>" +
    		"<span class='lihandle glyphicon' aria-hidden='true'></span>" +
    		"<span class='creator'></span>&nbsp;" +
    		"<span class='loc'></span>&nbsp;" +
    		"<span class='ctext'></span>" +
    		"</li>");
    
    if (this.creator == cm.getUser()) this.me.find(".creator").html("<em>You</em>");
    else this.me.find(".creator").text(this.creator);
    
    var handle = this.me.find(".lihandle");
    //change icon depending on privacy setting
    switch (this.ispublic)
    {
    case cm.PRIVACY_ONLYME:
    	handle.addClass("glyphicon-eye-close p_onlyme");
    	break;
    case cm.PRIVACY_PRIVATE:
    	handle.addClass("glyphicon-eye-open p_private");
    	break;
    case cm.PRIVACY_PUBLIC:
    	handle.addClass("glyphicon-globe p_public");
    	break;
    }
    
    this.me.find(".loc").text(this.chr+":"+this.start+".."+this.end);
	this.me.find(".ctext").text(this.ctext);
    
    this.deleteme = function() {
    	this.me.remove();
    };
    
    this.onclick = function(event) 
	{
    	event.stopPropagation();
    	cmdlg(this);
    	
    	if (this.isnew)
		{
    		this.me.removeClass("ne");
    		this.setnew(-1);
		}
	};
	this.me.on("click",this.onclick.bind(this));
}

/**
 * The structural node in the comment tree visualisation. The comment node can exist in 3 levels,
 * the root node - of which there can be 1, the set of track nodes - 1 for each custom data track,
 * and the chromosome node - for each chromosome with a registered comment. Chromosome nodes have Comment
 * objects as children while other nodes have CommentNode objects as children.
 */
function CommentNode(level, trackref, chr)
{
	this.indom = false;
	
	this.level = level;
	this.trackref = trackref;
	this.chr = chr;
	
	this.parent = null;
	this.children = {};
	this.numnewchildren = 0;
	
	//The node structure
	this.me = $("<section>" +
			    "<h4 class='ch'>" +
			    "	<span class='tex glyphicon glyphicon-menu-right' aria-hidden='true'></span>" +
			    "	<span class='cnt'></span>" +
			    "	<span class='ni'></span>" +
			    "</h4>" +
			    "<ul></ul>" +
			    "</section>");
	this.root = this.me.find('ul');
	this.handle = this.me.find('h4');
	this.expander = this.me.find('.tex');
	this.newindicator = this.me.find('.ni');
	this.deleteme = function() {
		this.me.remove();
	};

	this.expanded = false;
	this.content = this.me.find('.cnt');
	if (trackref == 0 && level == 0) 
	{
		this.expanded = true;
		this.expander.removeClass("glyphicon-menu-right");
		this.expander.addClass("glyphicon-menu-down");
		this.content.text("All Comments");
	}
	else if (level == 0)
	{
		this.content.text("TrackID: "+cm.tierIDMap[this.trackref]);
	}
	else
	{
		this.content.text("Chr"+this.chr+":");
		this.me.css("display","none");
	}
	
	/**
	 * Click event which expands or contracts the comment node depending on its current state
	 */
	this.onclick = function(event) 
	{
		event.stopPropagation();
		
		for (var child in this.children)
		{
			if (this.expanded)
			{
				this.children[child].me.hide();
			}
			else
			{
				this.children[child].me.show();
			}
		}
		this.expanded = !this.expanded;
		
		if (this.expanded == true)
		{
			this.expander.removeClass("glyphicon-menu-right");
			this.expander.addClass("glyphicon-menu-down");
		}
		else
		{
			this.expander.removeClass("glyphicon-menu-down");
			this.expander.addClass("glyphicon-menu-right");
		}
	};
	this.handle.on("click",this.onclick.bind(this));

	/**
	 * Function controlling the 'new children' notification, for when a new comment is made
	 * which is one of the children of the current node. Updates propagate from the bottom, up.
	 */
	this.newchildren = function(dif)
	{
		this.numnewchildren += dif;
		
		this.newindicator.html(this.numnewchildren);
		if (this.numnewchildren == 0) this.newindicator.hide();
		else this.newindicator.show(500);
		
		if (this.parent != undefined)
		{
			this.parent.newchildren(dif);
		}
	};
	
	/**
	 * Adds a new child to the current comment node
	 */
	this.addc = function (m,c)
	{
		this.children[m] = c;
		this.root.append(c.me);
	};
}

/**
 * The data descriptor object contains all the metadata properties of any given custom
 * data file. It also provides interaction via click events, allowing for various file
 * management options. Group-based file management and user file management, both utilise
 * the DD but perform different functions.
 */
function DataDesc(data, row, parenttable, option) 
{
    this.id 			= data[0];
    this.owner 			= data[1];
    this.trackname		= data[2];
    this.description	= data[3];
    this.path 			= data[4];
    this.format 		= data[5];
    this.ispublic 		= data[6];
    this.indexedby 		= data[7];
    this.isindex 		= data[8];
    this.remote 		= data[9];
    
    if (this.remote)
	{
    	this.meta = JSON.parse(this.format);
	}
    
    this.asArray = function()
    {
    	return [this.id, this.owner, 
    	        this.trackname, this.description, this.path, 
    	        (this.remote?JSON.stringify(this.meta):this.format), 
    	        this.ispublic, this.indexedby, 
    	        this.isindex, this.remote];
    };
    
    /**
     * 0 = file manager
     * 1 = public file manager
     * 2 = group access manager
     */
    this.foption = option;
    
    this.row = row;
    this.me = $(row);
    this.parent = parenttable;
    
    /**
     * Click event function defining how the file can be loaded as a track, as well as providing
     * means to manage the file, all via a dialog interface
     */
    this.def_onclick = function(e) 
	{
    	if (!e.metaKey && !e.ctrlKey && !e.shiftKey)
		{
    		//block dt selection unless shift/ctrl keys are pressed
    		e.preventDefault();
    		e.stopPropagation();
		}
    	else
		{
    		return;
		}

    	var selected = this.parent.rows({selected:true}).data();
    	var delete_offsets = this.parent.rows({selected:true})[0];
    	var datas = [];
    	var ids = [];
    	if (selected.length > 0)
		{
    		var this_selected = false;
    		
    		for (var i=0; i<selected.length; i++)
			{
    			datas.push(selected[i]);
    			ids.push(selected[i][0]);
    			
    			if (this.id == selected[i][0])
				{
    				this_selected = true;
				}
			}
    		
    		//only if this row is part of the selected group, get all selected rows
    		if (!this_selected)
			{
    			datas = [this.asArray()];
			}
		}
    	else
		{
    		datas.push(this.asArray());
    		ids.push(this.id);
		}
    	
    	//clear selected rows
    	this.parent.rows().deselect();
    	
    	var all_json = []; 
    	
    	for (var ind in datas)
		{
    		var meta = undefined;
    		if (datas[ind][9] == true)//this.remote)
    		{
    	    	meta = JSON.parse(datas[ind][5]);
    		}
    		
    		var json = {};
    		
    		//common to all track types
    		json["name"] = datas[ind][2]; //this.trackname;
    		json["desc"] = datas[ind][3]; //this.description;
    		
    		var fmt = undefined;
    		var pb_url = "//"+window.location.host+"/files/"+datas[ind][0];//this.id;
    		if (meta != undefined)
    		{
    			fmt = meta.data_format;
    			pb_url = datas[ind][4];//this.path;
    			if (meta.pennant)
    				json["pennant"] = meta.pennant;
    		}
    		else
			{
    			fmt = datas[ind][5];
			}
    		
    		switch (fmt.toLowerCase())
    		{
    		//same setup for bigwigs and bigbeds
    		case "bw":
    		case "bb":
    		case "bigwig":
    		case "bigbed":
    			json["bwgURI"] = pb_url;
    			json["collapseSuperGroups"] = true;
    			break;
    			
    		case "2bit":
    			json["twoBitURI"] = pb_url;
    			json["tier_type"] = "sequence"; 
    			break;
    			
    		case "bam":
    			json["bamURI"] = pb_url;
    			break;
    			
    		case "vcf":
    			json["uri"] = pb_url;
    			json["tier_type"] = "tabix";
    			json["payload"] = "vcf";
    			break;
    			
    		case "bed":
    			json["uri"] = pb_url;
    			json["tier_type"] = "tabix";
    			json["payload"] = "bed";
    			break;
    			
    		case "wig":
    			json["uri"] = pb_url;
    			json["tier_type"] = "memstore";
    			json["payload"] = "wig";
    			break;
    			
    		default:
    			isdata = false;
    			break;	
    		}
    		json["_ispublic"] = datas[ind][6];//this.ispublic;
    		json["_owner"] = datas[ind][1];//this.owner;
    		json["_id"] = datas[ind][0];//this.id;
    		
    		if (meta != undefined)
    		{
    			if (meta.feature_seq != "none")
    				json["tier_type"] = meta.feature_seq;
    		}
    		
    		all_json.push(json);
		}
    	
    	var fmdlg = $( "#file-manager-dlg" );
		fmdlg.find("#fmdlg_response").html("Select an option:");
		
		if (this.foption == 2)
		{
			//group file management
			$("#fm_standarddiv").addClass("hidden");
			$("#fm_groupdiv").removeClass("hidden");
		}
		else
		{
			$("#fm_standarddiv").removeClass("hidden");
			$("#fm_groupdiv").addClass("hidden");
		}
		
		var removefileaclbtn = fmdlg.find("#fm_removefileaclbtn");
		removefileaclbtn.off("click");
		removefileaclbtn.on("click", function(event) 
		{
			for (var id in ids)
			{
				cm.removeGroupFileACL(cm.active_group, ids[id]);
			}
			fmdlg.modal("hide");
		}.bind(this));
		
		var addtrackbtn = fmdlg.find("#fm_addtrackbtn");
		addtrackbtn.off("click");
		addtrackbtn.on("click", function(event) 
		{
			for (var js in all_json)
			{
				b.makeTier(all_json[js]);
			}
			b.markSelectedTiers();
			b.positionRuler();
			b.notifyTier();
			
			fmdlg.modal("hide");
		}.bind(this));
		
		var setonlyme = fmdlg.find("#fm_privacy_onlyme");
		setonlyme.off("click");
		setonlyme.on("click", function(event) 
		{
			for (var id in ids)
			{
				cm.togglePublicDataFile(ids[id],cm.PRIVACY_ONLYME);
			}
		}.bind(this));
		var setprivate = fmdlg.find("#fm_privacy_private");
		setprivate.off("click");
		setprivate.on("click", function(event) 
		{
			for (var id in ids)
			{
				cm.togglePublicDataFile(ids[id],cm.PRIVACY_PRIVATE);
			}
		}.bind(this));
		var setpublic = fmdlg.find("#fm_privacy_public");
		setpublic.off("click");
		setpublic.on("click", function(event) 
		{
			for (var id in ids)
			{
				cm.togglePublicDataFile(ids[id],cm.PRIVACY_PUBLIC);
			}
		}.bind(this));
		
		var alltoggles = $([]).add(setonlyme).add(setprivate).add(setpublic);
		alltoggles.removeClass("btn-primary");
		alltoggles.addClass("btn-default");
		switch (this.ispublic)
		{
		case cm.PRIVACY_ONLYME:
			setonlyme.addClass("btn-primary");
			setonlyme.removeClass("btn-default");
			break;
		case cm.PRIVACY_PRIVATE:
			setprivate.addClass("btn-primary");
			setprivate.removeClass("btn-default");
			break;
		case cm.PRIVACY_PUBLIC:
			setpublic.addClass("btn-primary");
			setpublic.removeClass("btn-default");
			break;
		}
		
		var delbtn = fmdlg.find("#fm_deletebtn");
		delbtn.off("click");
		delbtn.on("click", function(event) 
		{
			for (var id in ids)
			{
				cm.deleteDataFile(ids[id]);
			}
			if (ids.length == 1)
			{
				this.deleteme();
			}
			else
			{
				this.parent.row(delete_offsets).remove().draw();
			}
		}.bind(this));
		
		var sharetogroupbtn = fmdlg.find("#fm_sharetogroup");
		sharetogroupbtn.off("click");
		sharetogroupbtn.on("click", function(event) 
		{
			for (var id in ids)
			{
				cm.shareFileToGroup(ids[id], $('#fm_groupname').text());
			}
		}.bind(this));
		
		//can only do owner operations if all selected files are owned
		var ownall = true;
		for (var js in all_json)
		{
			if (cm.getUser() != all_json[js]["_owner"])
			{
				ownall = false;
			}
		}
		
		//don't show ineffectual options
		if (!ownall || this.foption == 1)
		{
			$("#fm_ownerops").hide();
			$("#fm_sharefile").hide();
		}
		else
		{
			$("#fm_ownerops").show();
			$("#fm_sharefile").show();
			
			//fill the select box
			var groups = fmdlg.find('#fm_grouplist');
			groups.html("");
			for (var i in cm.grouplist)
			{
				if ( (cm.grouplist[i]["access_level"]&2) != 2)
				{
					//give share option only for groups we have FILE_MANAGE privilege
					continue;
				}
				
				var li = $("<li><a href='#'>"+i+"</a></li>");
				li.on("click",function(){
					$('#fm_groupname').text($(this).find('a').text());
				});
				groups.append(li);
			}
		}
		
		fmdlg.modal('show');
		//find existing modal-backdrop
		var modalbd = $(".modal-backdrop:visible");
		modalbd[modalbd.length-1].style.zIndex = fmdlg.css("z-index")-1;
		
		//keep modal scrolling active
		fmdlg.on('hidden.bs.modal', function(event) {
			if (modalbd.length >= 1)
				$("body").addClass("modal-open");
		});
    };
    this.me.on("click",this.def_onclick.bind(this));
    
    //only the ispublic parameter can change
    this.update = function() {
    	var ddhtml = "";
    	
    	if (this.ispublic == cm.PRIVACY_PUBLIC) //PUBLIC
    		ddhtml = '<span class="glyphicon glyphicon-globe p_public" title="Everyone" aria-hidden="true"></span>';
    	else if (this.ispublic == cm.PRIVACY_PRIVATE) //PRIVATE
    		ddhtml = '<span class="glyphicon glyphicon-eye-open p_private" title="Yourself and Collaborators" aria-hidden="true"></span>';
		else if (this.ispublic == cm.PRIVACY_ONLYME) //ONLYME
    		ddhtml = '<span class="glyphicon glyphicon-eye-close p_onlyme" title="Only Yourself" aria-hidden="true"></span>';
    	
    	$('td', row).eq(6).html(ddhtml);
    };
    
    this.deleteme = function() {
    	this.parent.row(this.row).remove().draw();
    };
}

/**
 * The groupmember object allows for the performing of management tasks on group members.
 * Only users with the appropriate permissions may invoke these functions. Permissions are
 * all verified by the server upon receiving the request.
 */
function GroupMember(groupname, name, access_level, me, parenttable)
{
	this.groupname = groupname;
	this.name = name;
	this.access_level = access_level;
	
	/**
	 * Permissions are broken down into 4 types and stored as a set of flags in 1 integer.
	 * For example: the permission level 5 (in binary 0b0101) would indicate the user can
	 * manage users and read files i.e. 0b0001 | 0b0100 = 0b0101 = 5.
	 */
	this.flag_owner 		= (this.access_level&8)==8;
	this.flag_manage_users 	= (this.access_level&4)==4;
	this.flag_manage_files 	= (this.access_level&2)==2;
	this.flag_read 			= (this.access_level&1)==1;
	
	this.parent = parenttable;
	
	/**
	 * The user's permissions are visualised in html via a set of checkboxes
	 */
	this.me = $(me);
	this.row = me;
	this.acl = this.me.find(".gm_acl");
	
	this.setacl = function (acl)
	{
		this.access_level = acl;
		
		this.flag_owner 		= (this.access_level&8)==8;
		this.flag_manage_users 	= (this.access_level&4)==4;
		this.flag_manage_files 	= (this.access_level&2)==2;
		this.flag_read 			= (this.access_level&1)==1;
		
		this.me.find('.gm_p_r').prop('checked',this.flag_read);
		this.me.find('.gm_p_mf').prop('checked',this.flag_manage_files);
		this.me.find('.gm_p_mu').prop('checked',this.flag_manage_users);
		this.me.find('.gm_p_own').prop('checked',this.flag_owner);
		
		this.acl.text(acl);
	};
	this.setacl(this.access_level);

	this.invokeUpdateACL = function()
	{
		var p = 0;
		
		p += $(".mm_r")[0].checked?1:0;
		p += $(".mm_mf")[0].checked?2:0;
		p += $(".mm_mu")[0].checked?4:0;
		
		cm.updateGroupUserACL(this.groupname, this.name, p);
	};
	this.invokeDeleteUser = function() 
	{
		cm.removeUserFromGroup(this.groupname, this.name);
	};
	
	/**
	 * Click event for the groupmember object, creates a dialog box interface allowing managing
	 * users to perform administrative functions on them.
	 */
	this.onclick = function(event)
	{
		event.stopPropagation();
		
		if (this.access_level == 15 || (cm.grouplist[this.groupname]["access_level"]&4) != 4)
		{
			return;
		}
		
		var dlg = $("#group-user-dlg");
		dlg.find("#mm_acl").val(this.access_level);

		$(".mm_r")[0].checked = this.flag_read;
		$(".mm_mf")[0].checked = this.flag_manage_files;
		$(".mm_mu")[0].checked = this.flag_manage_users;
		
		$("#mm_setaclbtn").off("click");
		$("#mm_setaclbtn").on("click", function(event) 
		{
			dlg.modal("hide");
		});
		$("#mm_setaclbtn").on("click", this.invokeUpdateACL.bind(this));
		
		$("#mm_deleteuserbtn").off("click");
		$("#mm_deleteuserbtn").on("click", function(event) 
		{
			dlg.modal("hide");
		});
		$("#mm_deleteuserbtn").on("click", this.invokeDeleteUser.bind(this));
		
		dlg.modal('show');

		//find existing modal-backdrop
		var modalbd = $(".modal-backdrop:visible");
		modalbd[modalbd.length-1].style.zIndex = dlg.css("z-index")-1;
		
		//keep modal scrolling active
		dlg.on('hidden.bs.modal', function(event) {
			if (modalbd.length >= 1)
				$("body").addClass("modal-open");
		});
	};
	this.me.on("click",this.onclick.bind(this));
	
	this.deleteme = function() {
		this.parent.row(this.row).remove().draw();
    };
}

/**
 * The grouprow object is used to visualise the logged in user's group summary. Each object
 * represents 1 group which the user belongs to, and summarises the number of users and files
 * within it. It also shows the user's active permissions for that group.
 */
function GroupRow(groupname, num_members, num_files, access_level, me, parenttable)
{
	this.groupname = groupname;
	this.num_members = num_members;
	this.num_files = num_files;
	this.access_level = access_level;
	
	this.flag_owner 		= (this.access_level&8)==8;
	this.flag_manage_users 	= (this.access_level&4)==4;
	this.flag_manage_files 	= (this.access_level&2)==2;
	this.flag_read 			= (this.access_level&1)==1;
	
	this.parent = parenttable;
	
	this._row = me;
	this.row = $(me);

	this.updateFlags = function()
	{
		this.flag_owner 		= (this.access_level&8)==8;
		this.flag_manage_users 	= (this.access_level&4)==4;
		this.flag_manage_files 	= (this.access_level&2)==2;
		this.flag_read 			= (this.access_level&1)==1;
		
		this.row.find('.grow_p_r').prop('checked',this.flag_read);
		this.row.find('.grow_p_mf').prop('checked',this.flag_manage_files);
		this.row.find('.grow_p_mu').prop('checked',this.flag_manage_users);
		this.row.find('.grow_p_own').prop('checked',this.flag_owner);
	};
	//call it immediately
	this.updateFlags();
	
	/**
	 * Click event function which opens the detailed group management view.
	 * Depending on the user's permissions in that group, various management features may
	 * be unavailable to them.
	 */
	this.onclick = function(event)
	{
		//get data - we only need to call this once for each group
		//we will receive updates automatically afterwards
		if (cm.dd_group_data[this.groupname] == undefined)
		{
			cm.getGroupInfo(this.groupname);
		}
		
		//switch interface
		$("#groupoverviewdiv").addClass("hidden");
		$("#groupeditordiv").removeClass("hidden");
		$(".groupnameheader").text(this.groupname);
		
		if ((cm.grouplist[groupname]["access_level"]&8) == 8)
		{
			$("#deleteGroupBtn").removeClass("hidden");
			$("#leaveGroupBtn").addClass("hidden");
		}
		else 
		{
			$("#deleteGroupBtn").addClass("hidden");
			$("#leaveGroupBtn").removeClass("hidden");
		}
		
		this.hideMenus();
		
		//show group editor home button once we are browsing a group
		$(".editorbacknav").removeClass("hidden");
		
		//show entries from selected group
		cm.gfiletable.column(10).search(this.groupname).draw(false);
		cm.groupmembertable.column(1).search(this.groupname).draw(false);
		cm.active_group = this.groupname;
	};
	this.row.on("click",this.onclick.bind(this));
	
	this.hideMenus = function()
	{
		//selectively hide elements based on permissions
		if (!this.flag_manage_users)
			$("#ge_inu_tab").addClass("hidden");
		else $("#ge_inu_tab").removeClass("hidden");
		
		if (!this.flag_read && !this.flag_manage_files)
			$("#ge_rf_tab").addClass("hidden");
		else
			$("#ge_rf_tab").removeClass("hidden");
		
		if (!this.flag_manage_files)
			$("#ge_share_files").addClass("hidden");
		else $("#ge_share_files").removeClass("hidden");
		
		if (!this.flag_owner)
			$("#ge_gm_tab").addClass("hidden");
		else 
		{
			$("#ge_gm_tab").removeClass("hidden");
			
			$("#sendGroupMailBtn").off("click");
			$("#sendGroupMailBtn").on("click", function(event)
			{
				event.stopPropagation();
				var text = $("#groupemailtext").val();
				
				if (text.length == 0)
					return;
				
				cm.sendGroupMesssage(this.groupname, text);
				$("#groupemailtext").val("");
			}.bind(this));
		}
	};
	
	this.deleteme = function()
	{
		this.parent.row(this._row).remove().draw();
	};
	
	this.update = function()
	{
		this.row.find(".grow_gn").text(this.groupname);
		this.row.find(".grow_nm").text(this.num_members);
		this.row.find(".grow_nf").text((this.num_files==null?"0":this.num_files));
		this.row.find(".grow_al").text(this.access_level);
		
		this.updateFlags();
	};
}