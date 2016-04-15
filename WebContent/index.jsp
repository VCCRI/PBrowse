<%@ page language="java" contentType="text/html; charset=UTF-8"
		pageEncoding="UTF-8"%>
<html lang="en">

<head>

	<meta charset="utf-8">
	<meta http-equiv="X-UA-Compatible" content="IE=edge">
	<meta name="viewport" content="width=device-width, initial-scale=1">
	<meta name="description" content="">
	<meta name="author" content="PS">

	<title>PBrowse</title>

	<link href='https://fonts.googleapis.com/css?family=Roboto' rel='stylesheet' type='text/css'>

	<!-- Normalize CSS -->
	<link href="css/normalize.css" rel="stylesheet">

	<!-- DataTables CSS -->
	<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/s/bs/jq-2.1.4,dt-1.10.10,se-1.1.0/datatables.min.css"/>

	<!-- Bootstrap Core CSS -->
	<link href="css/bootstrap.min.css" rel="stylesheet">

	<!-- Custom CSS -->
	<link href="css/pbrowse.css" rel="stylesheet">

	<!-- jQuery -->
	<script src="js/jquery-2.1.4.min.js"></script>

	<!-- Bootstrap Core JavaScript -->
	<script src="js/bootstrap.min.js"></script>
	
	<!-- DataTables JS -->
	<script type="text/javascript" src="https://cdn.datatables.net/s/bs/dt-1.10.10,se-1.1.0/datatables.min.js"></script>

</head>

<body>

	<!-- Navigation -->
	<nav class="navbar navbar-inverse navbar-fixed-top" role="navigation">
		<div class="container">
			<!-- Brand and toggle get grouped for better mobile display -->
			<div class="navbar-header">
				<button type="button" class="navbar-toggle" data-toggle="collapse" data-target="#bs-example-navbar-collapse-1">
					<span class="sr-only">Toggle navigation</span>
					<span class="icon-bar"></span>
					<span class="icon-bar"></span>
					<span class="icon-bar"></span>
				</button>
				<a id="gv_menu" class="navbar-brand" href="#">PBrowse</a>
			</div>
			<!-- Collect the nav links, forms, and other content for toggling -->
			<div class="collapse navbar-collapse" id="bs-example-navbar-collapse-1">
				<ul class="nav navbar-nav">
					<li class="dropdown">
						<a href="#" class="dropdown-toggle" data-toggle="dropdown" role="button" aria-haspopup="true" aria-expanded="false">Genomes <span class="caret"></span></a>
						<ul class="dropdown-menu">
							<li class="dropdown-header">Generic</li>
							<li><a href="#" data-toggle="modal" data-target="#ggcDLG">Generic</a></li>
							<li role="separator" class="divider"></li>
							<li class="dropdown-header">H. sapiens (Human)</li>
							<li><a class="genome_selector" genome="hg38" href="#">GRCh38/hg38</a></li>
							<li><a class="genome_selector" genome="hg37" href="#">GRCh37/hg19</a></li>
							<li><a class="genome_selector" genome="ncbi36" href="#">NCBI36/hg18</a></li>
							<li role="separator" class="divider"></li>
							<li class="dropdown-header">M. musculus (Mouse)</li>
							<li><a class="genome_selector" genome="m38" href="#">GRCm38/mm10</a></li>
							<li><a class="genome_selector" genome="m37" href="#">NCBI37/mm9</a></li>
							<li role="separator" class="divider"></li>
							<li class="dropdown-header">D. rerio (Zebrafish)</li>
							<li><a class="genome_selector" genome="zv9" href="#">Zv9/danRer7</a></li>
							<li role="separator" class="divider"></li>
							<li class="dropdown-header">C. elegans (Worm)</li>
							<li><a class="genome_selector" genome="ws220" href="#">WS220/ce10</a></li>
						</ul>
			        </li>
			        <li class="dropdown">
						<a href="#" class="dropdown-toggle" data-toggle="dropdown" role="button" aria-haspopup="true" aria-expanded="false">Data <span class="caret"></span></a>
						<ul class="dropdown-menu">
							<li class="disabled"><a id="fu_menu" href="#" data-toggle="modal" data-target="#fileuploadDLG">File Upload</a></li>
							<li><a id="fm_menu" href="#" data-toggle="modal" data-target="#filemanagerDLG">File Management</a></li>
						</ul>
			        </li>
			        <li class="dropdown">
						<a href="#" class="dropdown-toggle" data-toggle="dropdown" role="button" aria-haspopup="true" aria-expanded="false">Collaboration <span class="caret"></span></a>
						<ul class="dropdown-menu">
							<li class="disabled"><a id="gm_menu" href="#" data-toggle="modal" data-target="#groupmanagerDLG">Groups</a></li>
							<li role="separator" class="divider"></li>
							<li class="disabled"><a id="cs_menu" href="#" data-toggle="modal" data-target="#newSessionDLG">Create New Session</a></li>
							<li class="disabled"><a id="js_menu" href="#" data-toggle="modal" data-target="#joinSessionDLG">Join Session</a></li>
							<li class="disabled"><a id="ps_menu" href="#" data-toggle="modal" data-target="#publicsessionviewDlg">Public Sessions</a></li>
						</ul>
			        </li>
			        <li class="dropdown">
						<a href="#" class="dropdown-toggle" data-toggle="dropdown" role="button" aria-haspopup="true" aria-expanded="false">Configs <span class="caret"></span></a>
						<ul class="dropdown-menu">
							<li class="disabled"><a id="cfgs_menu" href="#" data-toggle="modal" data-target="#configSaveDLG">Save Current Config</a></li>
							<li class="disabled"><a id="cfgl_menu" href="#" data-toggle="modal" data-target="#configLoadDLG">Load Existing Config</a></li>
						</ul>
			        </li>
			        <li><a id="forceRefresh" href="#" data-toggle="tooltip" data-placement="bottom" title="Force refresh of all loaded tiers."><span style="font-size: 20px;" class="glyphicon glyphicon-refresh" aria-hidden="true"></span></a></li>
				</ul>
				<ul class="nav navbar-nav navbar-right">
			        <li class="dropdown">
						<a id="profilemenu" href="#" class="hidden navbar-brand dropdown-toggle" data-toggle="dropdown" role="button" aria-haspopup="true" aria-expanded="false"><span class="loggedInUser"></span>&nbsp;<span class="caret"></span></a>
						<ul class="dropdown-menu">
							<li><a href="#" data-toggle="modal" data-target="#profilemanagerDLG">My Profile</a></li>
							<li role="separator" class="divider"></li>
							<li><a id="logoutButton" href="#" style="display:none">Logout</a></li>
						</ul>
			        </li>
			        <li><a id="loginButton" href="#" data-toggle="modal" data-target="#loginDLG">Login</a></li>
					<li><a id="registerButton" href="#" data-toggle="modal" data-target="#registerDLG">Register</a></li>
				</ul>
			</div>
			<!-- /.navbar-collapse -->
		</div>
		<!-- /.container -->
	</nav>

	<!-- Page Content -->
	<div class="container">

		<div class="row">

			<div class="col-md-8">
				<div id="genomePanel" class="">
					<div class="thumbnail" style="border: 0px;position: relative;">
						<div id="followerDivergeMask" class="modal-backdrop in" style="display:none;z-index: 6666;position: absolute;opacity: 0.15;"></div>
						<div class="img-responsive" id='svgHolder'></div>
					</div>
				</div>
			</div>

			<div class="col-md-4">
				<div id="sessionPanel" class="hidden">
					<div class="panel-group" id="session-accordion" role="tablist" aria-multiselectable="true">
					  <div class="panel panel-default">
					    <div class="panel-heading" role="tab" id="session-accordion-title">
					      <h4 class="panel-title">
					        <a role="button" data-toggle="collapse" data-parent="#session-accordion" href="#session-accordion-body" aria-expanded="true" aria-controls="collapseOne">
					          Session <span id="sessionLeaderIndicator"></span>
					        </a>
					      </h4>
					    </div>
					    <div id="session-accordion-body" class="panel-collapse collapse in" role="tabpanel" aria-labelledby="session-accordion-title">
					      <div class="panel-body">
							<div id="activeSessionDiv" style="display:none">
								<label id="sessionNameLabel"></label>
								<table style="width: 100%;">
									<tr>
										<td>SID:&nbsp;<span id="psessionID">(none)</span><span id="psessionCode">(none)</span></td>
										<td style="width: 250px;">
											<div class="btn-group btn-group-justified" role="group">
												<div class="btn-group" role="group">
													<button id="leaveSessionBtn" class="btn btn-warning">Leave Session</button>
												</div>
												<div class="btn-group" role="group">
													<button id="destroySessionBtn" class="btn btn-danger">Destroy Session</button>
												</div>
											</div>
										</td>
									</tr>
								</table>
								<label>User List:</label>
								<table class="table table-striped table-hover" style="word-break: break-all;">
									<tbody id="userListDiv">
										<!-- <tr><td>Username(Nickname)</td><td>gen-pos</td><td>leaderStatus</td></tr> -->
									</tbody>
								</table>
							</div>
							<div id="adminControllerDiv" style="display:none">
								<table style="width:100%">
									<tr>
										<td style="width:50%;padding-right: 7px;">
											<div class="form-group">
												<form id="admin_blacklistuser_form">
													<label>Blacklist User &nbsp;</label><span id="blacklistTooltip" class="glyphicon glyphicon-list" data-toggle="tooltip" data-placement="right" title="&lt;Empty&gt;"></span>
													<div id="sessionBlacklistUserError" style="display:none" class="hidden"></div>
													<div class="input-group">
										                <input id="sessionBlacklistUsername" type="text" class="form-control">
										                <span class="input-group-btn">
										                	<button id="sessionBlacklistUserDeny" class="btn btn-danger">Deny</button>
										                	<button id="sessionBlacklistUserAllow" class="btn btn-primary">Allow</button>
										                </span>
										            </div>
												</form>
											</div>
										</td>
										<td style="width:50%;padding-left: 7px;">
											<div class="form-group">
												<form id="admin_inviteuser_form">
													<label>Invite User</label>
													<div id="sessionInviteUserError" style="display:none" class="hidden"></div>
													<div class="input-group">
										                <input id="sessionInviteUsername" type="text" class="form-control">
										                <span class="input-group-btn">
										                	<button id="sessionInviteUserBtn" class="btn btn-default">Invite</button>
										                </span>
										            </div>
												</form>
											</div>
										</td>
									</tr>
								</table>
								
								<label class="checkbox-inline">
									<input type="checkbox" id="checkFollowerDiverge" checked value=""> Allow Divergence
								</label>
								<label class="checkbox-inline hidden">
									<input type="checkbox" id="checkFollowerSync" checked value=""> Follower Sync
								</label>
								<label class="checkbox-inline">
									<input type="checkbox" id="checkPrivateSession" checked value=""> Private Session
								</label>
							</div>
						  </div>
					      </div>
					    </div>
					</div>
					</div>
				
				<div id="sessionChatPanel" style="z-index: 6969" class="hidden">
					<div class="panel-group" id="chat-accordion" role="tablist" aria-multiselectable="true">
					  <div class="panel panel-default">
					    <div class="panel-heading" role="tab" id="headingOne">
					      <h4 class="panel-title">
					        <a class="collapsed" role="button" data-toggle="collapse" data-parent="#chat-accordion" href="#chat-accordion-body" aria-expanded="true" aria-controls="collapseOne">
					          Session Chat
					        </a>
					        <span class='ni'>0</span>
					      </h4>
					    </div>
					    <div id="chat-accordion-body" class="panel-collapse collapse" role="tabpanel" aria-labelledby="headingOne">
					      <div class="panel-body" style="padding-top: 14px;padding-bottom: 0px;">
							<div class="form-group">
								<div class="form-control" style="height:150px; overflow-y:auto" id="schatArea"></div>
							</div>
							<form id="sendMessageForm">
								<div class="form-group">
								    <label for="messageField">Enter a message:</label>
									<input type="text" class="form-control" id="messageField">
								</div>
								<input class="hidden" type="submit" tabindex="-1">
							</form>
					      </div>
					    </div>
					  </div>
					</div>
				</div>
			
				<div class="panel-group" id="comment-accordion" role="tablist" aria-multiselectable="true">
				  <div class="panel panel-default">
				    <div class="panel-heading" role="tab" id="headingOne">
				      <h4 class="panel-title">
				        <a role="button" data-toggle="collapse" data-parent="#comment-accordion" href="#comment-accordion-body" aria-expanded="true" aria-controls="collapseOne">
				          Comments
				        </a>
				      </h4>
				    </div>
				    <div id="comment-accordion-body" class="panel-collapse collapse in" role="tabpanel" aria-labelledby="headingOne">
				      <div class="panel-body">
<!-- 						<button id="newComment" disabled class="btn btn-default btn-block" style="width:100%;">New Comment</button> -->
<!-- 						<hr> -->
						<div id="commentlist">
						<!-- Automatically populated with nodes -->
						</div>
				      </div>
				    </div>
				  </div>
				</div>
				
				<!-- INCLUDE TESTS IF ?testing PARAM PROVIDED -->
				<% Boolean s = (Boolean) request.getAttribute("testmode"); %>
				<% if ( s != null ) { %>
				<%@ include file="/components/CollabTests.jsp" %>
				<% } %>
			</div>
		</div>
	</div>
	<!-- /.container -->

	<div class="pbrowse_footer">
		<div class="container">
			<hr>
			<!-- Footer -->
			<footer>
				<div class="row">
					<div class="col-lg-12">
						<p>Victor Chang Cardiac Research Institute 2016.</p>
					</div>
				</div>
			</footer>
		</div>
	</div>
	<!-- /.container -->

	<!-- Biodalliance Genome Browser core -->
	<script src="dalliance-compiled.js"></script>

	<!-- Biodalliance Genome Browser embed code -->
	<script src="js/GenomeBrowser.js"></script>

	<!-- Collabmanager -->
	<script src="js/HelperObjects.js"></script>

	<!-- Collabmanager -->
	<script src="js/CollabManager.js"></script>

	<!-- Dialog boxes -->
	<%@ include file="/components/Dialogs.jsp" %>

	<!-- Primary view JS -->
	<%@ include file="/components/IndexJS.jsp" %>
</body>

</html>
