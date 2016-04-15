<!--
UI related scripts
Function that works from non-dynamic (static) element
-->
<script>
	window.onbeforeunload = function() {
		//prevent refresh if the user is logged in with an active connection
		if (cm.getUser() != undefined && cm.sock.readyState == 1) {
			return "Are you sure you want to refresh? Session connectivity will be lost.";
		}
	};

	$(function() 
	{
		//registration dialog controls
		var emailRegex = /^[a-zA-Z0-9.!#$%&'*+\/=?^_`{|}~-]+@[a-zA-Z0-9](?:[a-zA-Z0-9-]{0,61}[a-zA-Z0-9])?(?:\.[a-zA-Z0-9](?:[a-zA-Z0-9-]{0,61}[a-zA-Z0-9])?)*$/; 

		//file upload
		var datafilename = undefined, datafileext = undefined, 
		indexfilename = undefined, indexfileext = undefined;

		function doLogin() {
			var vl = true;
			
			if ($("#l_username").val().length == 0 || $("#l_password").val().length == 0)
			{
				vl = false;
			}
			
			if (vl) {
				cm.doLogin($("#l_username").val(), $("#l_password").val());
			}
			return vl;
		}

		//registers a new user, if the input passes regex checks
		function addUser() {
			var valid = true;

			if (!(/^[a-z]([0-9a-z_\s])+$/i.test($("#username").val())))
			{
				$("#usernameerror").addClass("has-error");
				err($("#registerError"), "Username may consist of a-z, 0-9, underscores, spaces and must begin with a letter.", "danger");
				return;
			}
			else if ($("#username").val().length <=3 || $("#username").val().length > 12)
			{
				$("#usernameerror").addClass("has-error");
				err($("#registerError"), "Username must be between 4 and 12 characters.", "danger");
				return;
			}
			else $("#usernameerror").removeClass("has-error");
			
			if (!(/^[a-z]([0-9a-z_\s])+$/i.test($("#nickname").val())))
			{
				$("#nicknameerror").addClass("has-error");
				err($("#registerError"), "Nickname may consist of a-z, 0-9, underscores, spaces and must begin with a letter.", "danger");
				return;
			}
			else if ($("#nickname").val()<=3 || $("#nickname").val().length > 12)
			{
				$("#nicknameerror").addClass("has-error");
				err($("#registerError"), "Nickname must be between 4 and 12 characters.", "danger");
				return;
			}
			else $("#nicknameerror").removeClass("has-error");
			
			if (!(emailRegex.test($("#email").val())))
			{
				$("#emailerror").addClass("has-error");
				err($("#registerError"), "Invalid email, must follow: user@domain.com", "danger");
				return;
			}
			else $("#emailerror").removeClass("has-error");

			if ($("#password1").val().length < 5)
			{
				$("#passworderror").addClass("has-error");
				err($("#registerError"), "Password must be at least 5 characters", "danger");
				return;
			}
			else $("#passworderror").removeClass("has-error");
				
			if ($("#password1").val() != $("#password2").val())
			{
				$("#password2error").addClass("has-error");
				err($("#registerError"), "Passwords do not match", "danger");
				return;
			}
			else $("#password2error").removeClass("has-error");
			
			if (valid) {
				cm.doRegister($("#username").val(), $("#nickname").val(), $("#email").val(),
						$("#password1").val());
			}
			return valid;
		}

		//comment-form
		var cform = $("#comment-form").find("form").on("submit", function(event) {
			event.preventDefault();
		});
		var cdialog = $("#comment-form");
		cdialog.on('hide.bs.modal',function(){
			cform[0].reset();
		});
		cdialog.find('#cf_createcomment').on('click',function(e)
		{
			cm.makeComment($("#c_trackid").val(), $("#c_chr")
					.val(), $("#c_start").val(), $("#c_end")
					.val(), $("#c_comment").val(), $("#c_ispublic")
					.find("input:checked").val());

			cdialog.modal('hide');
		});
		//create session form
		var jsform = $('#newSessionDLG').find("form").on("submit",
			function(event) {
				event.preventDefault();
				cm.newSession();
		});
		//create session dialog opening event
		$('#newSessionDLG').on('show.bs.modal', function(event) {
			var modal = $(this);
			var newSession = modal.find("#newSessionBtn");
			jsform[0].reset();
			newSession.off("click");
			newSession.on("click", function() {
				cm.newSession();
			});
		});
		//file manager dialog opening event
		$('#filemanagerDLG').on('shown.bs.modal', function(event) {
			var modalbd = $(".modal-backdrop:visible");
			modalbd[modalbd.length-1].style.zIndex = $(this).css("z-index")-1;
			
			//keep modal scrolling active
			$(this).on('hidden.bs.modal', function(event) {
				if (modalbd.length > 1)
					$("body").addClass("modal-open");
			});
		});
		
		//join session form
		var nsform = $('#joinSessionDLG').find("form").on(
			"submit",
			function(event) {
				event.preventDefault();
				cm.joinSession($("#sessionIDText").val(), $(
						"#joinSessionCode").val());
			});
		//join session dialog opening event
		$('#joinSessionDLG').on(
			'show.bs.modal',
			function(event) {
				var modal = $(this);
				var joinSession = modal.find("#joinSessionBtn");
				nsform[0].reset();
				joinSession.off("click");
				joinSession.on("click", function() {
					cm.joinSession($("#sessionIDText").val(), $(
							"#joinSessionCode").val());
				});
			});
		//login form
		var lform = $('#loginDLG').find("form").on("submit", function(event) {
			event.preventDefault();
			doLogin();
		});
		//login dialog opening event
		$('#loginDLG').on('show.bs.modal', function(event) {
			var modal = $(this);
			var loginbtn = modal.find("#dologinbtn");
			var resetbtn = modal.find("#doresetbtn");
			lform[0].reset();

			loginbtn.off("click");
			loginbtn.on("click", function() {
				doLogin();
			});
			
			resetbtn.off("click");
			resetbtn.on("click", function() {
				cm.doPasswordReset();
			});
		});
		
		//reset fields in upload dialog
		$("#fileuploadDLG").on('show.bs.modal',function(){
			$("#fileUploadForm")[0].reset();
			$("#rm_fileUploadForm")[0].reset();
		});

		//register form
		var rform = $('#registerDLG').find("form").on("submit",
				function(event) {
					event.preventDefault();
					addUser();
				});
		//register dialog opening event
		$('#registerDLG').on('show.bs.modal', function(event) {
			var modal = $(this);
			var registerbtn = modal.find("#doregisterbtn");
			rform[0].reset();

			registerbtn.off("click");
			registerbtn.on("click", function() {
				addUser();
			});
		});


		$("#public-session-dlg").find("form").on("submit", function(event) {
			event.preventDefault();
		});

		$("#createGroupForm").on("submit", function(event) {
			event.preventDefault();
			$('#creategroupbutton').click();
		});

		$("#group-user-dlg").find("form").on("submit", function(event) {
			event.preventDefault();
		});

		var chatform = $("#sendMessageForm");
		chatform.on("submit", function(event) {
			event.preventDefault();
			var msg = $("#messageField").val();
			if (msg.length > 0 && msg.length < 256) {
				cm.sendSessionMsg(msg,false);
				chatform[0].reset();
			}
		});

		var regx = /^([0-9a-z_\s])+$/i;

		$('#creategroupbutton').click(function()
		{
			var gn = $("#groupname").val();
			if (regx.test(gn)) {
				cm.createGroup(gn);
				$("#createGroupError").hide();
			} else {
				err($("#createGroupError"),
						"The groupname can consist of alphanumeric characters only.","warning");
			}
			$("#createGroupForm")[0].reset();
		});

		//deleteGroupBtn
		$('#deleteGroupBtn').click(function() {
			var gn = cm.active_group;
			if (regx.test(gn)) {
				cm.deleteGroup(gn);
			}
		});

		//leaveGroupBtn
		$('#leaveGroupBtn').click(function() {
			var gn = cm.active_group;
			if (regx.test(gn)) {
				cm.leaveGroup(gn);
			}
		});

		//group navigation support
		var ge_all_tabs = $([]).add($("#groupmembers")).add($("#ge_invite_new_users_tab")).add($("#ge_group_mail_tab")).add($("#ge_read_files_tab"));
		$("#ge_invite_new_users").on("click",function(){
			ge_all_tabs.addClass("hidden");
			$("#ge_invite_new_users_tab").removeClass("hidden");
		});
		$("#ge_read_files").on("click",function(){
			ge_all_tabs.addClass("hidden");
			$("#ge_read_files_tab").removeClass("hidden");
		});
		$("#ge_group_mail").on("click",function(){
			ge_all_tabs.addClass("hidden");
			$("#ge_group_mail_tab").removeClass("hidden");
		});
		$(".groupbacknav").on("click", function(event) {
			event.stopPropagation();

			//revert tab view to group member list
			$('a[href="#groupmembers"]').tab('show');
			
			$("#groupoverviewdiv").removeClass("hidden");
			$("#groupeditordiv").addClass("hidden");
			$(".editorbacknav").addClass("hidden");
		});
		$(".editorbacknav").on("click", function(event) {
			event.stopPropagation();

			$("#groupoverviewdiv").addClass("hidden");
			$("#groupeditordiv").removeClass("hidden");
			
			ge_all_tabs.addClass("hidden");
			$("#groupmembers").removeClass("hidden");
		});

		$("#addNewUserBtn").on("click", function() {
			var p = 0;
			//calculate numerical permission value
			p += $("#nu_p_r")[0].checked ? 1 : 0; // file read
			p += $("#nu_p_f")[0].checked ? 2 : 0; // file manage
			p += $("#nu_p_u")[0].checked ? 4 : 0; // user manage
			cm.addUserToGroup($("#newUserUsername").val(), cm.active_group, p);
		});

		$("#logoutButton").on("click", function() {
			cm.doLogout();
		});

		//file upload menu
		$("#fu_menu").on("click", function(e) {
			if (this.parentNode.className == "disabled") e.stopPropagation();
		});
		//file manager menu
		$("#fm_menu").on("click", function(e) {
			if (this.parentNode.className == "disabled") e.stopPropagation();
		});
		//public session menu
		$("#ps_menu").on("click", function(e) {
			if (this.parentNode.className == "disabled") e.stopPropagation();
			else cm.getPublicSessions();
		});
		//group manager menu
		$("#gm_menu").on("click", function(e) {
			if (this.parentNode.className == "disabled") e.stopPropagation();
		});
		//create session menu
		$("#cs_menu").on("click", function(e) {
			if (this.parentNode.className == "disabled") e.stopPropagation();
		});
		//join session menu
		$("#js_menu").on("click", function(e) {
			if (this.parentNode.className == "disabled") e.stopPropagation();
		});
		//Save config menu
		$("#cfgs_menu").on("click", function(e) {
			if (this.parentNode.className == "disabled") e.stopPropagation();
		});
		//Load config menu
		$("#cfgl_menu").on("click", function(e) {
			if (this.parentNode.className == "disabled") e.stopPropagation();
		});

		function resetFormElement(e) {
			e.wrap('<form>').closest('form').get(0).reset();
			e.unwrap();
		}

		var require_index = false;

		//on change to the selected file for upload
		$(':file').change(function() {
			var file = this.files[0];
			
			var input = $(this).parents('.input-group').find(':text');
			input.val(file.name);
			
			var validexts = [ "bigbed", "bb", "bigwig", "bw", "2bit", "bam", "vcf",
					"bed", "wig", "bai", "tbi" ];
			var indexable = [ "bam", "vcf", "bed" ];
			var validdata = [ "bigbed", "bb", "bigwig", "bw", "2bit", "bam", "vcf",
					"bed", "wig" ];
			var validindex = [ "bai", "tbi" ];
			var exts = file.name.toLowerCase().split(".");
			var ext = "";
			for (var i = exts.length - 1; i >= 0; i--) {
				for (var j = 0; j < validexts.length; j++) {
					if (exts[i] == validexts[j]) {
						ext = validexts[j];
						break;
					}
				}
				if (ext != "")
					break;
			}
			if ($(this).hasClass("dataUploadFile")) {
				if (file == undefined) {
					datafilename = null;
					return;
				}
				datafilename = file.name;

				datafileext = ext;

				if (validdata.indexOf(ext) == -1) {
					err($("#fileuploadError"), "Not a valid genomic data file.", "danger");
					$("#inputFileError").addClass("has-error");
					input.val("");
					resetFormElement($(this));
				} else if (indexable.indexOf(ext) != -1) {
					$(".indexUploadFile").show();
					require_index = true;
				} else {
					require_index = false;					
					$(".indexUploadFile").hide();
				}
			} else if ($(this).hasClass("indexUploadFile")) {
				if (file == undefined) {
					datafilename = null;
					return;
				}

				indexfilename = null;
				indexfileext = ext;
				if (validindex.indexOf(ext) == -1) {
					err($("#fileuploadError"), "Not a valid index file!", "danger");
					$(".indexUploadFile").addClass("has-error");
					input.val("");
					resetFormElement($(this));
					return;
				} else if (datafileext == "bam") {
					if (indexfileext != "bai") {
						err($("#fileuploadError"), "Expected a .bai file!", "danger");
						input.val("");
						resetFormElement($(this));
						return;
					}
				} else if (datafileext == "vcf" || datafileext == "bed") {
					if (indexfileext != "tbi") {
						err($("#fileuploadError"), "Expected a .tbi file!", "danger");
						input.val("");
						resetFormElement($(this));
						return;
					}
				}
				$(".indexUploadFile").removeClass("has-error");
				$("#inputFileError").removeClass("has-error");
				
				indexfilename = file.name;
			}
		});

		var progressbar = $("#progressbar");
		var pbarlabel = $("#progressbarlabel");
		progressbar.css("width",'0%');
		
		var formData = null;
		$('#uploadFileButton').click(function() {
			if (cm.getUser() == undefined) {
				return;
			}
			
			var valid = true;

			if (!(/^([0-9a-z_])+$/i.test($("#studyid").val())))
			{
				$("#studyiderror").addClass("has-error");
				err($("#fileuploadError"), "Study ID may consist of a-z, 0-9, and underscores.", "danger");
				return;
			}
			else $("#studyiderror").removeClass("has-error");
			
			if (!(/^([0-9a-z_\s])+$/i.test($("#trackname").val())))
			{
				$("#tracknameerror").addClass("has-error");
				err($("#fileuploadError"), "Trackname may consist of a-z, 0-9, underscores, and spaces.", "danger");
				return;
			}
			else $("#tracknameerror").removeClass("has-error");
			
			if (!(/^([0-9a-z_\s])*$/i.test($("#description").val())))
			{
				$("#descriptionerror").addClass("has-error");
				err($("#fileuploadError"), "Description may consist of a-z, 0-9, underscores, and spaces.", "danger");
				return;
			}
			else $("#descriptionerror").removeClass("has-error");
			
			valid = valid && (datafilename != null);

			if (require_index && indexfilename == null) {
				$("#indexUploadFile").addClass("has-error");
				valid = false;
				err($("#fileuploadError"), "Index file not selected!", "warning");
			}
			else $("#indexUploadFile").removeClass("has-error");

			if (!valid)
				return;

			formData = new FormData($('#fileUploadForm')[0]);
			$.ajax({
				url : 'uploadData', //Server script to process data
				type : 'POST',
				xhr : function() { // Custom XMLHttpRequest
					var myXhr = $.ajaxSettings.xhr();
					if (myXhr.upload) { // Check if upload property exists
						myXhr.upload
							.addEventListener(
								'progress',
								progressHandlingFunction,
								false); // For handling the progress of the upload
					}
					return myXhr;
				},
				//Ajax events
				success : completeHandler,
				error : errorHandler,
				// Form data
				data : formData,
				//Options to tell jQuery not to process data or worry about content-type.
				cache : false,
				contentType : false,
				processData : false
			});
			
			progressbar.addClass("active");

			$(this).hide();
		});

		function completeHandler(e) {
			var msg = JSON.parse(e);
			if (msg.status == "failed") {
				err($("#fileuploadError"), msg.error, "danger");
			} else {
				err($("#fileuploadError"), "File upload successful", "success");
				$("#fileUploadForm")[0].reset();
			}
			$("#uploadFileButton").show();
			progressbar.removeClass("active");
		}
		function errorHandler(e) {
			err($("#fileuploadError"), "File upload failed.", "danger");
		}

		function progressHandlingFunction(e) {
			if (e.lengthComputable) {
				var width = (((e.loaded+1)/(e.total+1))*100) | 0;
				progressbar.css("width",width+'%');
				pbarlabel.text(width+'%');
			}
		}

		//remote data-source registration ------------------------------
		$("a[href=#localDataTab]").on("click", function(){
			$("#rm_uploadFileButton").addClass("hidden");
			$("#uploadFileButton").removeClass("hidden");
		});
		$("a[href=#remoteDataTab]").on("click", function(){
			$("#rm_uploadFileButton").removeClass("hidden");
			$("#uploadFileButton").addClass("hidden");
		});
		
		$("#rm_dataurl").on("input",function() {
			var url = $(this).val();
			console.log("change");
			
			if (url.length == 0)
			{
				$("#rm_trackname").val("");
				$("#rm_studyid").val("");
			}
			
			if ( !(/^(http|https):\/\/[a-z0-9\.\_\-\%\/]+$/i.test(url)) )
			{
				$("#rm_dataurl_error").addClass("has-error");
			}
			else 
			{
				$("#rm_dataurl_error").removeClass("has-error");
				var p = url.split(".");
				var ext = p[p.length-1];
				
				p = url.split("/");
				var filename = p[p.length-1];
				
				var fn = filename.split(".");
				var combined_fn = "";
				for (var i=0; i<fn.length-1; i++)
				{
					combined_fn += fn[i];
					if (i != fn.length-2 && fn.length-2>0)
						combined_fn+=".";
				}
				$("#rm_trackname").val(combined_fn);
				
				var studyid = url.split("/")[2];
				$("#rm_studyid").val(studyid);
				
				$("#rm_feature_seq").parent().addClass("hidden");
				
				//if it matches known formats, we can autofill some fields
				switch (ext.toLowerCase())
				{
				case "bam":
					$("#rm_data_format").find("[value='bam']")[0].checked = true;
					$("#rm_feature_seq").find("[value=none]")[0].checked = true;
					break;

				case "bb":
				case "bigbed":
					$("#rm_data_format").find("[value='bigbed']")[0].checked = true;
					$("#rm_feature_seq").find("[value=none]")[0].checked = true;
					break;
				case "bw":
				case "bigwig":
					$("#rm_data_format").find("[value='bigwig']")[0].checked = true;
					$("#rm_feature_seq").find("[value=none]")[0].checked = true;
					break;
					
				case "vcf":
					$("#rm_feature_seq").parent().removeClass("hidden");
					$("#rm_data_format").find("[value='vcf']")[0].checked = true;
					$("#rm_feature_seq").find("[value=tabix]")[0].checked = true;
					break;
					
				case "bed":
					$("#rm_feature_seq").parent().removeClass("hidden");
					$("#rm_data_format").find("[value='bed']")[0].checked = true;
					$("#rm_feature_seq").find("[value=tabix]")[0].checked = true;
					break;
					
				case "wig":
					$("#rm_feature_seq").parent().removeClass("hidden");
					$("#rm_data_format").find("[value='wig']")[0].checked = true;
					$("#rm_feature_seq").find("[value=memstore]")[0].checked = true;
					break;
					
				case "2b":
				case "2bit":
					$("#rm_feature_seq").parent().removeClass("hidden");
					$("#rm_data_format").find("[value='2bit']")[0].checked = true;
					$("#rm_feature_seq").find("[value=sequence]")[0].checked = true;
					break;
					
				//assume this is a DAS tier
				default:
					$("#rm_feature_seq").parent().removeClass("hidden");
					$("#rm_data_format").find("[value='das']")[0].checked = true;
					$("#rm_feature_seq").find("[value=sequence]")[0].checked = true;
					break;
				}
				
			}
		});
		
		$("#rm_uploadFileButton").on("click",function() {
			var url = $("#rm_dataurl").val();
			if ( !(/^(http|https):\/\/[a-z0-9\.\_\-\%\/]+$/i.test(url)) )
			{
				err($("#remoteuploadError"),"Invalid URL provided.","danger");
				$("#rm_dataurl_error").addClass("has-error");
				return;
			} else $("#rm_dataurl_error").removeClass("has-error");
			
			var studyid = $("#rm_studyid").val();
			if (!(/^[a-z0-9\.]+$/i.test(studyid)))
			{
				err($("#remoteuploadError"),"Study id must consist of letters, numbers and fullstops.","danger");
				$("#rm_studyid").parent().addClass("has-error");
				return;
			} else $("#rm_studyid").parent().removeClass("has-error");
			
			var trackname = $("#rm_trackname").val();
			if (!(/^[a-z0-9\.]+$/i.test(trackname)))
			{
				err($("#remoteuploadError"),"Trackname must consist of letters, numbers and fullstops.","danger");
				$("#rm_trackname").parent().addClass("has-error");
				return;
			} else $("#rm_trackname").parent().removeClass("has-error");
			
			cm.registerRemoteSource($("#rm_dataurl").val(), 
					$("#rm_data_format").find("[name=data_format]:checked").val(),
					$("#rm_feature_seq").find("[name=feature_seq]:checked").val(), 
					$("#rm_studyid").val(), $("#rm_trackname").val(),
					$("#rm_genomic_class").find("[name=genometag]:checked").val(),
					$("#rm_description").val(), 
					$("#rm_toggleispublic").find("[name=ispublic]:checked").val()
			);
		});
		//--------------------------------------------------------------
		
		
		//genome-viewer menu: Collab session manager
		$("#leaveSessionBtn").on("click", function() {
			cm.leaveSession();
		});

		$("#destroySessionBtn").on("click", function() {
			cm.endSession();
		});

		//admin panel checkboxes
		$("#checkFollowerSync").change(function() {
			cm._enableFollowerSync(this.checked);
		});
		$("#checkFollowerDiverge").change(function() {
			cm._allowDivergance(this.checked);
		});
		$("#checkPrivateSession").change(function() {
			cm._setPrivateSession(this.checked);
		});

		$(".genome_selector").on("click", function() {
			cm.changeGenome($(this).attr("genome"), function(){
				b.notifyTier();
			});
		});
		
		$("#forceRefresh").on("click",function(){
			//reload current browser - fixes cached loading errors workaround
			cm.changeGenome(cm.active_genome);
		});

		$("#profilemanagerDLG").find("form").each(function() {
			$(this).on("submit", function(event) {
				event.preventDefault();
				$(this)[0].reset();
			});
		});

		$("#changePasswordBtn").on(
			"click",
			function() {
				if ($("#pm_newpassword").val().length < 5)
				{
					err($("#changePasswordError"),"New password must be at least 5 characters!","danger");
					return;
				}
				if ($("#pm_oldpassword").val().length == 0)
				{
					err($("#changePasswordError"),"Fields cannot be left empty.","info");
					return;
				}
				if ($("#pm_newpassword").val() != $("#pm_newpassword2").val())
				{
					err($("#changePasswordError"),"New passwords do not match!","danger");
					return;
				}
				
				$("#changePasswordError").addClass("hidden");
				cm.changeUserPassword($("#pm_oldpassword").val(), 
						$("#pm_newpassword").val(), $("#pm_newpassword2").val());
			}
		);

		$("#changeEmailBtn").on("click", function() {
			if (!emailRegex.test($("#pm_newemail").val()))
			{
				err($("#changeEmailError"),"Email address should be of the form: name@domain.com","danger");
				return;
			}
			$("#changeEmailError").addClass("hidden");
			cm.changeUserEmail($("#pm_newemail").val());
		});

		$("#changeNicknameBtn").on("click", function() {
			if ($("#pm_newnickname").val().length < 3 || $("#pm_newnickname").val().length > 12)
			{
				err($("#changeNicknameError"),"Nickname must be least 3-12 characters","danger");
				return;
			}
			$("#changeNicknameError").addClass("hidden");
			cm.changeUserNickname($("#pm_newnickname").val());
		});
		
		//enchanced leader option controls
		$("#sessionBlacklistUserAllow").on("click", function(e) {
			e.preventDefault();e.stopPropagation();
			var name = $("#sessionBlacklistUsername").val();
			if (name.length > 0)
				cm.sessionBlacklistUpdate(name,"allow");
		});
		$("#sessionBlacklistUserDeny").on("click", function(e) {
			e.preventDefault();e.stopPropagation();
			var name = $("#sessionBlacklistUsername").val();
			if (name.length > 0)
				cm.sessionBlacklistUpdate(name,"deny");
		});
		$("#sessionInviteUserBtn").on("click", function(e) {
			e.preventDefault();e.stopPropagation();
			var name = $("#sessionInviteUsername").val();
			if (name.length > 0)
				cm.sessionInviteUser(name);
		});
		//invite user submit form on enter behaviour
		$("#admin_inviteuser_form").on("submit", function(e){
			e.preventDefault();e.stopPropagation();
			$("#sessionInviteUserBtn").click();
		});
		//blacklist user default behaviour on enter
		$("#admin_blacklistuser_form").on("submit", function(e){
			e.preventDefault();e.stopPropagation();
			$("#sessionBlacklistUserDeny").click();
		});
		
		//indicator control for session chat
		$("#sessionChatPanel").find("a[href=#chat-accordion-body]").on("click",function(){
			var ind = $("#sessionChatPanel").find(".ni");
			ind.text("0");
			ind.hide(100);
			
			//modify body padding when box is open
			if ($("#sessionChatPanel").find(".collapsed").length > 0)
				$("body").css("padding-bottom","353px");
			else
				$("body").css("padding-bottom","93px");
		});
		
		//save config button
		$("#cfg_save_btn").on("click", function(e){
			e.preventDefault();
			
			if (cm.getUser() == undefined)
			{
				return;
			}
			
			var name = $("#cfg_save_name").val();
			var desc = $("#cfg_save_description").val();
			
			var genomename = cm.active_genome;
			if (cm.active_genome == "generic")
			{
				genomename = genomes["generic"].coordSystem.speciesName;
			}
			
			//add the entry to the saved config list
			cm.savedConfigTable.row.add(
				[name,desc,genomename,new Date().toLocaleString(),JSON.stringify(cm.getTiers())]
			).draw();
			
			//push update
			cm.updateSavedConfig();
		});
		
		//GENOME CREATOR CONTROLS ~~~~~~~~~~~~~~~~~~~~~~~~~~~
		$("#ggc_create_btn").on("click",function(e){
			
			var name = $("#ggc_genome_name").val();
			var seqfile = $("#ggc_seq_file").val();
			var genefile = $("#ggc_gene_file").val();
			
			if (seqfile.length == 0)
			{
				$("#ggc_seq_file").parent().addClass("has-error");
				err($("#ggc_form_error"),"Sequence file is required.","danger");
				return false;
			}
			if (name.length == 0)
			{
				err($("#ggc_form_error"),"Genome name is required.","danger");
				$("#ggc_genome_name").parent().addClass("has-error");
				return false;
			}
			if ($("#ggcDLG").find(".has-error").length > 0)
			{
				err($("#ggc_form_error"),"Fix all errors before proceeding.","danger");
				return false;
			}
			
			var gen = genomes["generic"];
			
			var f = undefined;
			if (cm.dd_private_data[seqfile] != undefined)
			{
				f = cm.dd_private_data[seqfile];
			}
			else if (cm.dd_public_data[seqfile] != undefined)
			{
				f = cm.dd_public_data[seqfile];
			}
			else
			{
				$("#ggc_seq_file").parent().addClass("has-error");
				err($("#ggc_form_error"),"Invalid sequence file ID.","danger");
				return false;
			}
			var path = "//"+window.location.host+"/files/"+seqfile;
			if (f.remote)
			{
				path = f.path;
			}
			
			var g = undefined;
			if (genefile.length > 0)
			{
				if (cm.dd_private_data[genefile] != undefined)
				{
					g = cm.dd_private_data[genefile];
				}
				else if (cm.dd_public_data[genefile] != undefined)
				{
					g = cm.dd_public_data[genefile];
				}
				else
				{
					$("#ggc_gene_file").parent().addClass("has-error");
					err($("#ggc_form_error"),"Invalid gene file ID.","danger");
					return false;
				}
				
				path = "//"+window.location.host+"/files/"+genefile;
				if (f.remote)
				{
					path = g.path;
				}
				
				gen.sources.push({
					name: g.trackname,
					desc: g.description,
					bwgURI: path,
					stylesheet_uri: '//www.biodalliance.org/stylesheets/gencode2.xml',
					trixURI: path+".ix",
					collapseSuperGroups: true,
				});
			}

			gen.coordSystem.speciesName = name;
			//add the file as a source
			gen.sources.push({
				name: f.trackname,
				desc: f.description,
				twoBitURI: path,
				tier_type: 'sequence',
				_ispublic: f.ispublic,
				_owner: f.owner,
				_id: seqfile,
			});
			
			//load the new genome as the generic template
			cm.changeGenome("generic", function() {
				$('#ggcDLG').modal('hide');
				b.notifyTier();
			});
		});
		$("#ggc_seq_file").on("input",function() {
			if (!(/^[0-9]+$/.test($(this).val())))
			{
				$(this).parent().addClass("has-error");
			} 
			else $(this).parent().removeClass("has-error");
		});
		$("#ggc_genome_name").on("input",function() {
			if (!(/^[0-9a-z\._]+$/i.test($(this).val())))
			{
				$(this).parent().addClass("has-error");
			} 
			else $(this).parent().removeClass("has-error");
		});
		$("#ggc_gene_file").on("input",function() 
		{
			if (!(/^[0-9]*$/.test($(this).val())))
			{
				$(this).parent().addClass("has-error");
			} 
			else $(this).parent().removeClass("has-error");
		});
		$('#ggcDLG').on('shown.bs.modal', function(event) {
			var dlg = $(this);
			//find existing modal-backdrop
			var modalbd = $(".modal-backdrop:visible");
			modalbd[modalbd.length-1].style.zIndex = dlg.css("z-index")-1;
			
			//keep modal scrolling active
			dlg.on('hidden.bs.modal', function(event) {
				if (modalbd.length >= 1)
					$("body").addClass("modal-open");
			});
			
			dlg.find('form')[0].reset();
		});
		$("#ggcDLG").find("form").on("submit", function(e) {
			e.preventDefault();e.stopPropagation();
			$("#ggc_create_btn").click();
		});
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		$("#configSaveDLG").on('shown.bs.modal', function(e) {
			var dlg = $(this);
			//find existing modal-backdrop
			var modalbd = $(".modal-backdrop:visible");
			modalbd[modalbd.length-1].style.zIndex = dlg.css("z-index")-1;
			
			//keep modal scrolling active
			dlg.on('hidden.bs.modal', function(event) {
				if (modalbd.length >= 1)
					$("body").addClass("modal-open");
			});
			
			dlg.find('form')[0].reset();
		});
		$("#configSaveDLG").find("form").on("submit",function(e){
			e.stopPropagation();e.preventDefault();
		});
		
		//disable anchor link behaviour
		$("body").find("a").each(function(){
			if (this.href.indexOf("#") != -1)
			{
				$(this).on("click",function(ev){ ev.preventDefault(); });
			}
		});
		
		//capture f5 event, changing it to cause a tier refresh through our code
		$(document.body).on("keydown", this, function (event) {
		    if (event.keyCode == 116) {
		    	cm.changeGenome(cm.active_genome);
		        event.stopPropagation();
		        event.preventDefault();
		    }
		});
		
		//enable tooltips
		$('[data-toggle="tooltip"]').tooltip();
	});

	function showAlert(title, msg) {
		$("#alertDlgLabel").text(title);
		$("#alertDlgText").html(msg);
		var aldlg = $("#alertDlg");
		aldlg.modal({
			  keyboard: false,
			  show: true,
			  backdrop: 'static'
		});
		
		var mbd = $(".modal-backdrop:visible");
		mbd[mbd.length-1].style.zIndex = aldlg.css("z-index")-1;
		aldlg.on('hidden.bs.modal',function(){
			if (mbd.length > 1)
				$("body").addClass("modal-open");
		});
	}

	//for viewing and removing comments
	function cmdlg(c) {
		//clear it
		var com = $("#cmdlg_comment");
		com.text(c.ctext);
		$("#cmdlg_position").val(c.chr + ":" + c.start + ".." + c.end);
		$("#cmdlg_author").val(c.creator);

		var _cmdlg = $("#cmdlg");

		var delbtn = _cmdlg.find('#cmdlg_del_btn');
		delbtn.addClass("hidden"); //by default
		//allow the creator to delete, only
		if (cm.getUser() == c.creator) {
			delbtn.removeClass("hidden");
			delbtn.off("click");
			delbtn.on("click",function(){
				cm.deleteComment(c);
				cm._removeComment(c);
				_cmdlg.modal("hide");
			});
		}

		var gotobtn = _cmdlg.find('#cmdlg_goto_btn');
		gotobtn.off("click");
		gotobtn.on("click",function(){
			b.setLocation(c.chr, c.start, c.end);
			_cmdlg.modal("hide");
		});

		_cmdlg.modal('show');
	}

	function toggleLoginInterfaces(bool) {
		var bars = $([]).add("#fu_menu").add("#gm_menu")
		.add("#cs_menu").add("#js_menu").add("#ps_menu")
		.add("#cfgs_menu").add("#cfgl_menu");
		if (!bool) bars.each(function(){$(this).parent().addClass("disabled");});
		else bars.each(function(){$(this).parent().removeClass("disabled");});
	}

	function err(ele, msg, status) {
		ele.removeClass("hidden");
		ele.text(msg);
		ele.addClass("alert alert-"+status);
		ele.show(500);
		setTimeout(function() {
			ele.hide(500,function(){ ele.removeClass("alert alert-"+status); });
		}, 5000);
	}
</script>
