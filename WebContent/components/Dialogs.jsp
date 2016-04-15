<!-- File Management Dialog BS - Public and Private -->
<div class="modal" id="filemanagerDLG" style="z-index: 9100 !important;" tabindex="-1" role="dialog" aria-labelledby="filemanagerDLGLabel">
	<div class="modal-dialog modal-lg" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="filemanagerDLGLabel">Manage Files</h4>
			</div>
			<div class="modal-body">
				<ul class="nav nav-tabs" role="tablist">
					<li role="presentation" class="active"><a href="#privateDataTab" aria-controls="privateDataTab" role="tab" data-toggle="tab">My Dataset</a></li>
					<li role="presentation"><a href="#publicDataTab" aria-controls="publicDataTab" role="tab" data-toggle="tab">Public Dataset</a></li>
				</ul>
				<!-- Tab panes -->
				<div class="tab-content" style="margin-top: 12px;">
					<div role="tabpanel" class="tab-pane active" id="privateDataTab">
						<p style="text-align: justify;">Below are the files you have previously uploaded to PBrowse. Use the interface here to add any of them as tracks
						in the genome-viewer interface. If you are part of a collaborative session, all tracks listed here as excluding <em>only-me</em> tracks, will be
						shared with all participants for the duration of the session.</p>
						<p>Shift+click or Ctrl+Click to select multiple entries at once.</p>
						<div id="fileManagerDiv">
							<table id="fileTable" style="width:100%;" class="table table-striped table-bordered table-hover">
								<thead>
									<tr><th>ID</th><th>Owner</th><th>Track-name</th>
									<th>Description</th>
									<th>Path</th><th>Format</th>
									<th>Public</th><th>Index ID</th>
									<th>_isIndex</th><th>_isRemote</th></tr>
								</thead>
								<tbody>
								
								</tbody>
							</table>
						</div>
					</div>
					<div role="tabpanel" class="tab-pane" id="publicDataTab">
						<p style="text-align: justify;">Below are the files which have been uploaded for public use by any users of pbrowse. Use the interface here to add any of them as tracks
						in the genome-viewer interface. If you are part of a collaborative session, any tracks enabled here will be
						shared with all participants for the duration of the session.</p>
						<p>Shift+click or Ctrl+Click to select multiple entries at once.</p>
						<div id="publicfileManagerDiv">
							<table id="publicfileTable" style="width:100%;" class="table table-striped table-bordered table-hover">
								<thead>
									<tr><th>ID</th><th>Owner</th><th>Track-name</th>
									<th>Description</th>
									<th>Path</th><th>Format</th>
									<th>Public</th><th>Index ID</th>
									<th>_isIndex</th><th>_isRemote</th></tr>
								</thead>
								<tbody>
								
								</tbody>
							</table>
						</div>
					</div>
				</div>
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
			</div>
		</div>
	</div>
</div>

<!-- File Upload Dialog BS -->
<div class="modal" id="fileuploadDLG" tabindex="-1" role="dialog" aria-labelledby="fileuploadDLGLabel">
	<div class="modal-dialog modal-lg" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="fileuploadDLGLabel">Upload Files</h4>
			</div>
			<div class="modal-body">
				<ul class="nav nav-tabs" role="tablist">
					<li role="presentation" class="active"><a href="#localDataTab" aria-controls="localDataTab" role="tab" data-toggle="tab">Upload Data</a></li>
					<li role="presentation"><a href="#remoteDataTab" aria-controls="remoteDataTab" role="tab" data-toggle="tab">Remote Data</a></li>
				</ul>
				<!-- Tab panes -->
				<div class="tab-content" style="margin-top: 12px;">
					<div role="tabpanel" class="tab-pane active" id="localDataTab">
						<p style="text-align: justify;">Select a file you wish to upload for viewing with the genome browser - currently supports <em>wig, 
						bed<sup>1</sup>, bigbed, bigwig, vcf, 2bit and bam</em> formats. For indexes, the <em>bai and tbi</em> formats are allowed,
						depending on input data type. If your data file requires an index, you will be prompted to upload it at the same time.</p>
						<p style="text-align: justify; font-size:12px;">* : Bed file needs to be saved in bgzipped format (.bed.bz) and require a tabix index.</p>
						<p style="text-align: justify;">Control who can view your uploads, or delete uploaded files by visiting the file manager.</p>
						<div id="fileuploadError" style="display:none" class="hidden"></div>
						<form enctype="multipart/form-data" id="fileUploadForm" role="form">
							
							<div class="form-group">
								<label for="data_file">Data file:</label>
								<span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="The data file to upload. You will be prompted for an index if required."></span>
								<div id="inputFileError" class="input-group">
					                <span class="input-group-btn">
					                    <span class="btn btn-default btn-file">
					                        Browse&hellip; <input class="dataUploadFile" name="data_file" type="file">
					                    </span>
					                </span>
					                <input type="text" style="background-color: #FFF;" class="form-control" readonly>
					            </div>
							</div>
							<div class="form-group">
								<div class="indexUploadFile" style="display:none;">
									<label for="index_file">Index file:</label>
						            <span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="The index file associated with the selected data file. The correct index format for the given file must be chosen."></span>
								</div>
								<div style="display:none;" class="input-group indexUploadFile">
					                <span class="input-group-btn">
					                    <span class="btn btn-default btn-file">
					                        Browse&hellip; <input name="index_file" class="indexUploadFile" type="file">
					                    </span>
					                </span>
					                <input type="text" style="background-color: #FFF;" class="form-control" readonly>
					            </div>
							</div>
							<div id="studyiderror" class="form-group">
								<label for="studyid">Study identifier:</label>
								<span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="Data categorisation, for grouping of related files."></span>
								<input type="text" name="studyid" id="studyid" placeholder="myfirststudy" class="form-control">
							</div>
							<div id="tracknameerror" class="form-group">
								<label for="trackname">Track-name:</label>
								<span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="The name of the track when shown in the browser."></span>
								<input type="text" name="trackname" id="trackname" placeholder="track01" class="form-control">
							</div>
							<div class="form-group">
								<label>Species:</label>
								<span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="Description of the species. User compliance is optional."></span>
								<div id="genomic_class" style="margin-bottom: 8px;">
									<label class="radio-inline"><input type="radio" value="#Human" name="genometag" checked>Human</label>
									<label class="radio-inline"><input type="radio" value="#Mouse" name="genometag">Mouse</label>
									<label class="radio-inline"><input type="radio" value="#Zebrafish" name="genometag">Zebrafish</label>
									<label class="radio-inline"><input type="radio" value="#Worm" name="genometag">Worm</label>
<!-- 									<label class="radio-inline"><input type="radio" value="#Aspergillus" name="genometag">A. Nidulans</label> -->
								</div>
							</div>
							<div id="descriptionerror" class="form-group">
								<label for="description">Description:</label>
								<span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="Extra data related to the track. Information here can be searched and filtered."></span>
								<input type="text" name="description" id="description" placeholder="This is the track description" class="form-control">
							</div>
							<div class="form-group">
								<label>File Visibility:</label>
								<span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="Affects who can see the data, can be changed later on."></span>
								<div id="toggleispublic" style="margin-bottom: 8px;">
									<label class="radio-inline"><input type="radio" name="ispublic" value=1 checked>Only-Me</label>
									<label class="radio-inline"><input type="radio" name="ispublic" value=2 >Collaborators</label>
									<label class="radio-inline"><input type="radio" name="ispublic" value=4 >Public</label>
								</div>
							</div>
							<div class="progress">
								<div id="progressbar" class="progress-bar progress-bar-striped active" role="progressbar" style="width:0%">
									<span id="progressbarlabel">0%</span>
								</div>
							</div>
						</form>
					</div> <!-- END LOCAL DATA TAB -->
					<div role="tabpanel" class="tab-pane" id="remoteDataTab">
						<p style="text-align: justify;">If your data is already hosted somewhere else, register it with PBrowse by entering
						the required information below. It will be searchable with all your other files.</p>
						<p style="text-align: justify;">For data types requiring an index, ensure it can be reached at the same url with appropriate extension.</p>
						<div id="remoteuploadError" style="display:none" class="hidden"></div>
						<form id="rm_fileUploadForm" role="form">
							<div id="rm_dataurl_error" class="form-group">
								<label for="description">URL:</label>
								<span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="The url to the remotely hosted file. Must be accessible through http(s)."></span>
								<input type="text" id="rm_dataurl" name="dataurl" placeholder="http://server/myfile.sort.bam" class="form-control">
							</div>
							<div class="form-group">
								<label>Data-Format:</label>
								<span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="The type of data, according to its extension. Defaults to DAS if no extension."></span>
								<div id="rm_data_format">
									<label class="radio-inline"><input type="radio" value="bam" name="data_format" checked>BAM</label>
									<label class="radio-inline"><input type="radio" value="bigbed" name="data_format">BigBed</label>
									<label class="radio-inline"><input type="radio" value="bigwig" name="data_format">BigWig</label>
									<label class="radio-inline"><input type="radio" value="vcf" name="data_format">VCF</label>
									<label class="radio-inline"><input type="radio" value="bed" name="data_format">Bed</label>
									<label class="radio-inline"><input type="radio" value="wig" name="data_format">Wig</label>
									<label class="radio-inline"><input type="radio" value="2bit" name="data_format">2Bit</label>
									<label class="radio-inline"><input type="radio" value="das" name="data_format">DAS</label>
								</div>
							</div>
							<div class="form-group">
								<label>Data-Type:</label>
								<span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="The tier type for loading in Dalliance. Leave as-is if unsure."></span>
								<div id="rm_feature_seq">
									<label class="radio-inline"><input type="radio" value="none" name="feature_seq" checked>None</label>
									<label class="radio-inline"><input type="radio" value="sequence" name="feature_seq">Sequence</label>
									<label class="radio-inline"><input type="radio" value="tabix" name="feature_seq">Tabix</label>
									<label class="radio-inline"><input type="radio" value="memstore" name="feature_seq">Memstore</label>
								</div>
							</div>
							<div class="form-group">
								<label for="rm_studyid">Study identifier:</label>
								<span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="Data categorisation, for grouping of related files."></span>
								<input type="text" name="rm_studyid" id="rm_studyid" placeholder="myfirststudy" class="form-control">
							</div>
							<div class="form-group">
								<label for="rm_trackname">Track-name:</label>
								<span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="The name of the track when shown in the browser."></span>
								<input type="text" name="rm_trackname" id="rm_trackname" placeholder="track01" class="form-control">
							</div>
							<div class="form-group">
								<label>Species:</label>
								<span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="Description of the species. User compliance is optional."></span>
								<div id="rm_genomic_class" style="margin-bottom: 8px;">
									<label class="radio-inline"><input type="radio" value="#Human" name="genometag" checked>Human</label>
									<label class="radio-inline"><input type="radio" value="#Mouse" name="genometag">Mouse</label>
									<label class="radio-inline"><input type="radio" value="#Zebrafish" name="genometag">Zebrafish</label>
									<label class="radio-inline"><input type="radio" value="#Worm" name="genometag">Worm</label>
<!-- 									<label class="radio-inline"><input type="radio" value="#Aspergillus" name="genometag">A. Nidulans</label> -->
								</div>
							</div>
							<div id="rm_description_error" class="form-group">
								<label for="rm_description">Description:</label>
								<span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="Extra data related to the track. Information here can be searched and filtered."></span>
								<input type="text" name="rm_description" id="rm_description" placeholder="This is the track description" class="form-control">
							</div>
							<div class="form-group">
								<label>File Visibility:</label>
								<span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="Affects who can see the data, can be changed later on."></span>
								<div id="rm_toggleispublic" style="margin-bottom: 8px;">
									<label class="radio-inline"><input type="radio" name="ispublic" value=1 checked>Only-Me</label>
									<label class="radio-inline"><input type="radio" name="ispublic" value=2 >Collaborators</label>
									<label class="radio-inline"><input type="radio" name="ispublic" value=4 >Public</label>
								</div>
							</div>
						</form>
					</div> <!-- END REMOTE DATA TAB -->
				</div>
			
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
				<button id="uploadFileButton" type="button" class="btn btn-primary">Upload</button>
				<button id="rm_uploadFileButton" type="button" class="btn btn-primary hidden">Upload</button>
			</div>
		</div>
	</div>
</div>

<!-- Register Dialog BS -->
<div class="modal" id="registerDLG" tabindex="-1" role="dialog" aria-labelledby="registerDLGLabel">
	<div class="modal-dialog modal-sm" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="registerDLGLabel">Register</h4>
			</div>
			<div class="modal-body">
				<form id="registerForm">
					<div id="registerError" style="display:none" class="hidden"></div>
					<div id="usernameerror" class="form-group">
					  <label for="username">Username:</label>
					  <input type="text" class="form-control" id="username" placeholder="Jane Smith">
					</div>
					<div id="nicknameerror" class="form-group">
					  <label for="nickname">Nickname:</label>
					  <input type="text" class="form-control" id="nickname" placeholder="JaneS">
					</div>
					<div id="emailerror" class="form-group">
					  <label for="email">Email Address:</label>
					  <input type="email" class="form-control" id="email" placeholder="me@domain.com">
					</div>
					<div id="passworderror" class="form-group">
					  <label for="l_username">Password</label>
					  <input type="password" class="form-control" id="password1">
					</div>
					<div id="password2error" class="form-group">
					  <label for="l_password">Repeat Password</label>
					  <input type="password" class="form-control" id="password2">
					</div>
					<input type="submit" tabindex="-1" class="hidden">
				</form>
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
				<button id="doregisterbtn" type="button" class="btn btn-primary">Register</button>
			</div>
		</div>
	</div>
</div>

<!-- Login Dialog BS -->
<div class="modal" id="loginDLG" tabindex="-1" role="dialog" aria-labelledby="loginDLGLabel">
	<div class="modal-dialog modal-sm" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="loginDLGLabel">Login</h4>
			</div>
			<div class="modal-body">
				<form id="loginform">
					<div id="loginError" style="display:none" class="hidden"></div>
					<div class="form-group">
					  <label for="l_username">Username</label>
					  <input type="text" class="form-control" id="l_username" placeholder="Username">
					</div>
					<div class="form-group">
					  <label for="l_password">Password</label>
					  <input type="password" class="form-control" id="l_password">
					</div>
					<input type="submit" tabindex="-1" class="hidden">
				</form>
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
				<button id="doresetbtn" type="button" class="btn btn-warning">Reset Password</button>
				<button id="dologinbtn" type="button" class="btn btn-primary">Login</button>
			</div>
		</div>
	</div>
</div>

<!-- New Session Dialog BS -->
<div class="modal" id="newSessionDLG" tabindex="-1" role="dialog" aria-labelledby="nsdlg_Label">
	<div class="modal-dialog modal-sm" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="nsdlg_Label">New Session</h4>
			</div>
			<div class="modal-body">
				<form id="newsessionform">
					<div class="form-group">
					  <label for="newSessionName">Session Name:</label>
					  <input type="text" class="form-control" id="newSessionName" placeholder="MySession" autocomplete="off">
					</div>
					<div class="form-group">
					  <label for="newSessionCode">Password:</label>
					  <input type="password" class="form-control" id="newSessionCode" autocomplete="off">
					</div>
					<div class="checkbox" style="margin-top: 0px;">
						<label><input type="checkbox" value="" checked id="newSessionPrivate">Make private</label>
					</div>	
					<input type="submit" tabindex="-1" class="hidden">
				</form>
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
				<button id="newSessionBtn" type="button" class="btn btn-primary">Create Session</button>
			</div>
		</div>
	</div>
</div>

<!-- Join Session Dialog BS -->
<div class="modal" id="joinSessionDLG" tabindex="-1" role="dialog" aria-labelledby="jsdlg_Label">
	<div class="modal-dialog modal-sm" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="jsdlg_Label">Join Session</h4>
			</div>
			<div class="modal-body">
				<form id="joinsessionform">
					<div id="joinsessionError" style="display:none" class="hidden"></div>
					<div class="form-group">
					  <label for="sessionIDText">Session Name:</label>
					  <input type="text" class="form-control" id="sessionIDText" placeholder="ABCDEF" autocomplete="off">
					</div>
					<div class="form-group">
					  <label for="joinSessionCode">Password:</label>
					  <input type="password" class="form-control" id="joinSessionCode" autocomplete="off">
					</div>
					<input type="submit" tabindex="-1" class="hidden">
				</form>
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
				<button id="joinSessionBtn" type="button" class="btn btn-primary">Join Session</button>
			</div>
		</div>
	</div>
</div>

<!-- File Manager Dialog BS -->
<div class="modal" id="file-manager-dlg" style="z-index: 10001 !important;" tabindex="-1" role="dialog" aria-labelledby="fmdlg-Label">
	<div class="modal-dialog modal-fm" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="fmdlg-Label">File Options</h4>
			</div>
			<div class="modal-body">
				<div style="display:none" id="fmdlg_response">Select an option:</div>
				<div class="form-group">
					<button id="fm_addtrackbtn" class="btn btn-default btn-block">Add As Track</button>
				</div>
				<div id="fm_standarddiv">
					<div id="fm_ownerops">
						<div class="form-group">
							<div class="btn-group btn-group-justified" role="group">
								<div class="btn-group" role="group">
									<button type="button" id="fm_privacy_onlyme" class="btn btn-default">Only-Me</button>
								</div>
								<div class="btn-group" role="group">
									<button type="button" id="fm_privacy_private" style="padding-left: 2px; padding-right: 2px;" class="btn btn-default">Collaborators</button>
								</div>
								<div class="btn-group" role="group">
									<button type="button" id="fm_privacy_public" class="btn btn-default">Public</button>
								</div>
							</div>
						</div>
						<div class="form-group">
							<button id="fm_deletebtn" class="btn btn-danger btn-block">Delete File</button>
						</div>
					</div>
					<div id="fm_sharefile">
						<div class="btn-group" role="group" style="width: 100%;">
							<div class="btn-group" role="group" style="width: 75%;">
								<button type="button" class="btn btn-default dropdown-toggle" data-toggle="dropdown" style="width: 100%;" aria-haspopup="true" aria-expanded="false">
									<span id="fm_groupname">Select Group:</span>
									<span class="caret" style="position: absolute;right: 10px;top: 15px;"></span>
								</button>
								<ul id="fm_grouplist" class="dropdown-menu" style="width: 100%;text-align: center;">
									<li><a href="#">None</a></li>
								</ul>
							</div>
							<div class="btn-group" role="group" style="width: 25%;">
								<button id="fm_sharetogroup" type="button" class="btn btn-default" style="width: 100%;">Share</button>
							</div>
						</div>
					</div>
				</div>
				<div id="fm_groupdiv" style="margin-top: 6px;" class="hidden">
					<div class="form-group">
						<button id="fm_removefileaclbtn" class="btn btn-danger btn-block">Remove File Access</button>
					</div>
				</div>
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
			</div>
		</div>
	</div>
</div>

<!-- Session User Manager Dialog BS -->
<div class="modal" id="sum-dlg" style="z-index: 10001 !important;" tabindex="-1" role="dialog" aria-labelledby="fmdlg-Label">
	<div class="modal-dialog modal-sm" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="fmdlg-Label"><span id="sum_title"></span></h4>
			</div>
			<div class="modal-body">
				<div id="sum_leaderopts" class="hidden">
					<div class="form-group">
						<button id="sum_setleader" class="btn btn-primary btn-block">Set Leader</button>
					</div>
					<div class="form-group">
						<button id="sum_kickuser" class="btn btn-warning btn-block">Kick User</button>
					</div>
					<div class="form-group">
						<button id="sum_blacklistuser" class="btn btn-danger btn-block">Blacklist User</button>
					</div>
				</div>
				<div id="sum_followeropts" class="hidden">
					<p>Only the leader of the session is permitted to manage other collaborators.
					Select an option below to proceed.</p>
				</div>
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
				<button id="sum_follow_user" type="button" class="btn btn-default">Follow</button>
				<button id="sum_gotopos" type="button" class="btn btn-primary">Goto Position</button>
			</div>
		</div>
	</div>
</div>

<!-- Group Member Manager Dialog BS -->
<div class="modal" id="group-user-dlg" style="z-index: 10001 !important;" tabindex="-1" role="dialog" aria-labelledby="gudlg-Label">
	<div class="modal-dialog modal-sm" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="gudlg-Label">Manage Member</h4>
			</div>
			<div class="modal-body">
				<div class="form-group">
					<label>Select a combination of the following:</label>
					<div class='checkbox'>
					<label><input class='mm_r' type='checkbox' value=''>Read Files</label>
					</div>
					<div class='checkbox'>
					<label><input class='mm_mf' type='checkbox' value=''>Manage Files</label>
					</div>
					<div class='checkbox'>
					<label><input class='mm_mu' type='checkbox' value=''>Manage Users</label>
					</div>
					<button id="mm_setaclbtn" class="btn btn-default btn-block">Set ACL</button>
					<button id="mm_deleteuserbtn" class="btn btn-danger btn-block">Remove User</button>
				</div>
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
			</div>
		</div>
	</div>
</div>

<!-- Public Session Viewer Dialog BS -->
<div class="modal" id="publicsessionviewDlg" tabindex="-1" role="dialog" aria-labelledby="psvDlgLabel">
	<div class="modal-dialog modal-lg" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="psvDlgLabel">Public Sessions</h4>
			</div>
			<div class="modal-body">
				<p>Below are all the currently running, publicly visible sessions. Sessions without an entrycode may
				be joined via the collaborative menu panel in the genome-viewer menu.</p>
				<div id="publicSessionManagerDiv">
					<table style="width:100%;" class="table table-striped table-bordered table-hover" id="publicSessionTable">
						<thead>
							<tr><th>SID</th><th>Session Name</th><th>Number of Participants</th><th>Requires Entrycode</th></tr>
						</thead>
						<tbody>
						
						</tbody>
					</table>
				</div>
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
			</div>
		</div>
	</div>
</div>

<div class="modal" id="public-session-dlg" tabindex="-1" role="dialog" aria-labelledby="psDlgLabel">
	<div class="modal-dialog modal-sm" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="psDlgLabel">Public Session</h4>
			</div>
			<div class="modal-body">
				<div class="form-group">
					<form>
						<div class="form-group">
							<label for="ps_sid">Session ID:</label>
							<input type="text" name="ps_sid" id="ps_sid" class="form-control">
						</div>
						<div class="form-group" id="ps_passcode_row">
							<label for="ps_sid">Passcode:</label>
							<input type="text" name="ps_passcode" id="ps_passcode" class="form-control">
						</div>
					</form>
				</div>
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
				<button id="ps_joinbtn" class="btn btn-primary">Join Session</button>
			</div>
		</div>
	</div>
</div>

<!-- Alert Dialog BS -->
<div class="modal" id="alertDlg" style="z-index: 10010 !important;" tabindex="-1" role="dialog" aria-labelledby="alertDlgLabel">
	<div class="modal-dialog modal-sm" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="alertDlgLabel">Alert</h4>
			</div>
			<div class="modal-body">
				<table>
					<tr>
						<td style="width: 30px;"><span class="glyphicon glyphicon-alert" aria-hidden="true"></span></td>
						<td><span id="alertDlgText"></span></td>
					</tr>
				</table>
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
			</div>
		</div>
	</div>
</div>

<!-- Group Manager Dialog BS -->
<div class="modal" id="groupmanagerDLG" style="z-index: 9000 !important;" tabindex="-1" role="dialog" aria-labelledby="gmDLGLabel">
	<div class="modal-dialog modal-lg" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="gmDLGLabel">Group Management</h4>
			</div>
			<div class="modal-body">
				<div id="groupoverviewdiv">
					<ol class="breadcrumb">
						<li class="active">Groups</li>
					</ol>
					<p>Manage all the groups you have admin privileges with, join groups, create new groups, leave them, share files with them.</p>
					<label>Your Groups:</label>
					<table id="groupRoleTable" style="width:100%;" class="table table-striped table-bordered table-hover">
						<thead>
							<tr><th>Groupname</th><th>Number of Members</th>
							<th>Number of Files</th><th>R</th><th>F</th><th>U</th><th>O</th></tr>
						</thead>
						<tbody>
						
						</tbody>
					</table>
					<div id="createGroupError" style="display:none" class="hidden"></div>					
					<form id="createGroupForm">
						<div class="form-group">
							<label for="groupname">Create a new group:</label>
							<input type="text" style="margin-bottom: 4px;" name="groupname" id="groupname" placeholder="Groupname" class="form-control">
							<button id="creategroupbutton" class="btn btn-default btn-block">Create</button>
						</div>
					</form>
				</div>
				<div id="groupeditordiv" class="hidden">
					<ol class="breadcrumb">
						<li><a class="groupbacknav" href="#">Groups</a></li>
						<li class="active"><span class="groupnameheader"></span></li>
					</ol>
					<ul class="nav nav-tabs" role="tablist">
						<li role="presentation" class="active"><a href="#groupmembers" aria-controls="groupmembers" role="tab" data-toggle="tab">Group Members</a></li>
						<li role="presentation"><a id="ge_inu_tab" href="#ge_invite_new_users_tab" aria-controls="ge_invite_new_users_tab" role="tab" data-toggle="tab">Invite New Users</a></li>
						<li role="presentation"><a id="ge_rf_tab" href="#ge_read_files_tab" aria-controls="ge_read_files_tab" role="tab" data-toggle="tab">Shared Files</a></li>
						<li role="presentation"><a id="ge_gm_tab" href="#ge_group_mail_tab" aria-controls="ge_group_mail_tab" role="tab" data-toggle="tab">Send Group Email</a></li>
					</ul>
					<!-- Tab panes -->
					<div class="tab-content" style="margin-top: 12px;">
						<div role="tabpanel" class="tab-pane active" id="groupmembers">
							<label>Group Members:</label>
							<table id="groupMemberTable" style="width:100%;" class="table table-striped table-bordered table-hover">
								<thead>
									<tr><th>Name</th><th>_GroupName</th><th>Access Level</th><th>R</th><th>F</th><th>U</th><th>O</th></tr>
								</thead>
								<tbody>
								
								</tbody>
							</table>
							<div id="ge_share_files">
								<br>
								<label>Share a file with this group:</label>
								<p>To share a file with this group go to your <a href="#" data-toggle="modal" data-target="#filemanagerDLG">file manager</a> and select the
								file, then select <em>share with group</em>. Then choose "<span class="groupnameheader"></span>" from the list. You can remove access via the table 
								above or change the required access.</p>
							</div>
						</div>
						<div role="tabpanel" class="tab-pane" id="ge_invite_new_users_tab">
							<div id="addNewUserError" style="display:none" class="hidden"></div>
							<div class="form-group">
								<label for="newUserUsername">Username:</label>
								<input type="text" name="newUserUsername" id="newUserUsername" placeholder="username" class="form-control">
							</div>
							<div class="form-group">
								<label>Permissions:</label>
								<label class="checkbox-inline"><input type="checkbox" id="nu_p_r" value=""> Read Files</label>
								<label class="checkbox-inline"><input type="checkbox" id="nu_p_f" value=""> Manage Files</label>
								<label class="checkbox-inline"><input type="checkbox" id="nu_p_u" value=""> Manage Users</label>
							</div>
							<button id="addNewUserBtn" style="margin-top: 5px;margin-bottom: 20px;" class="btn btn-default btn-block">Add User</button>
						</div>
						<div role="tabpanel" class="tab-pane" id="ge_read_files_tab">
							<label>Shared Files:</label>
							<p>Shift+click or Ctrl+Click to select multiple entries at once.</p>
							<table style="width:100%;" class="table table-striped table-bordered table-hover" id="groupFileTable">
								<thead>
									<tr><th>ID</th><th>Owner</th>
									<th>Track-name</th><th>Description</th>
									<th>Path</th><th>Format</th><th>Is Public</th>
									<th>Indexed By</th><th>IsIndex</th><th>_isremote</th><th>_Groupname</th></tr>
								</thead>
								<tbody>
								
								</tbody>
							</table>
						</div>
						<div role="tabpanel" class="tab-pane" id="ge_group_mail_tab">
							<div id="groupemailError" style="display:none" class="hidden"></div>
							<div class="form-group">
								<label for="groupemailtext">Your message:</label>
								<textarea class="form-control" rows="5" id="groupemailtext"></textarea>
							</div>
							<button id="sendGroupMailBtn" class="btn btn-default btn-block">Send</button>
						</div>
						<div role="tabpanel" class="tab-pane" id="publicDataTab">
						</div>
					</div>
				</div>
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
				<button id="deleteGroupBtn" class="hidden btn btn-danger">Delete Group</button>
				<button id="leaveGroupBtn" class="hidden btn btn-danger">Leave Group</button>
			</div>
		</div>
	</div>
</div>

<!-- Profile Manager Dialog BS -->
<div class="modal" id="profilemanagerDLG" tabindex="-1" role="dialog" aria-labelledby="pmDLGLabel">
	<div class="modal-dialog modal-lg" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="pmDLGLabel">My Profile</h4>
			</div>
			<div class="modal-body">
				<p>Modify your user profile below:</p>
				<form id="changePasswordForm" role="form">
					<label>Change password:</label>
					<div id="changePasswordError" style="display:none" class="hidden"></div>
					<table style="width:100%">
						<tr>
							<td style="width: 120px; text-align: right; padding-right: 8px;">Old Password:</td>
							<td><input type="password" name="pm_oldpassword" id="pm_oldpassword" class="form-control"></td>
						</tr>
						<tr>
							<td style="width: 120px; text-align: right; padding-right: 8px;">New Password:</td>
							<td><input type="password" name="pm_newpassword" id="pm_newpassword" class="form-control"></td>
						</tr>
						<tr>
							<td style="width: 120px; text-align: right; padding-right: 8px;">Repeat Password:</td>
							<td><input type="password" name="pm_newpassword2" id="pm_newpassword2" class="form-control"></td>
						</tr>
						<tr>
							<td></td>
							<td><button id="changePasswordBtn" class="btn btn-default btn-block">Submit</button></td>
						</tr>
					</table>
				</form>
				<br>
				<form id="changeEmailForm" role="form">						
					<label>Change Associated Email:</label>
					<div id="changeEmailError" style="display:none" class="hidden"></div>
					<table style="width:100%">
						<tr>
							<td style="width: 120px; text-align: right; padding-right: 8px;">New Email:</td>
							<td><input type="text" name="pm_newemail" id="pm_newemail" placeholder="newemail@domain.com" class="form-control"></td>
						</tr>
						<tr>
							<td></td>
							<td><button id="changeEmailBtn" class="btn btn-default btn-block">Submit</button></td>
						</tr>
					</table>
				</form>	
				<br>					
				<form id="changeNicknameForm" role="form">
					<label>Change Nickname:</label>
					<div id="changeNicknameError" style="display:none" class="hidden"></div>
					<table style="width:100%">
						<tr>
							<td style="width: 120px; text-align: right; padding-right: 8px;">New Email:</td>
							<td><input type="text" name="pm_newnickname" id="pm_newnickname" placeholder="MyNickname" class="form-control"></td>
						</tr>
						<tr>
							<td></td>
							<td><button id="changeNicknameBtn" class="btn btn-default btn-block">Submit</button></td>
						</tr>
					</table>
				</form>
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
			</div>
		</div>
	</div>
</div>

<!-- Comment Form Dialog BS -->
<div class="modal" id="comment-form" tabindex="-1" role="dialog" aria-labelledby="cfDLGLabel">
	<div class="modal-dialog" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="cfDLGLabel">New Comment</h4>
			</div>
			<div class="modal-body">
				<form>
				<fieldset>
					<input type="hidden" id="c_trackid">
					<div class="form-group">
						<div class="form-inline">
							<div class="input-group input-group-sm">
								<div class="input-group-addon">ID</div>
								<input readonly type="text" style="background-color: #FFF;" class="form-control" id="c_trackname">
							</div>
							<div class="input-group input-group-sm" style="float:right">
								<div class="input-group-addon">chr</div>
								<input readonly type="text" style="background-color: #FFF;width: 35px;text-align: center;" class="form-control" id="c_chr">
								<div class="input-group-addon">:</div>
								<input readonly type="text" style="background-color: #FFF;width: 80px;text-align: center;" class="form-control" id="c_start">
								<div class="input-group-addon">-</div>
								<input readonly type="text" style="background-color: #FFF;width: 80px;text-align: center;" class="form-control" id="c_end">
							</div>
						</div>
					</div>
					<div class="form-group">
						<label>Visibility:</label>
						<div id="c_ispublic" style="margin-bottom: 8px;">
							<label class="radio-inline"><input name="c_ispublic_radio" type="radio" value=1 checked>Only-Me</label>
							<label class="radio-inline"><input name="c_ispublic_radio" type="radio" value=2 >Collaborators</label>
							<label class="radio-inline"><input name="c_ispublic_radio" type="radio" value=4 >Public</label>
						</div>
					</div>
					<div class="form-group">
						<label for="c_comment">Comment:</label>
						<textarea id="c_comment" class="form-control" rows="3" placeholder="Why was this region commented?"></textarea>
					</div>
				</fieldset>
				</form>
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
				<button id="cf_createcomment" type="button" class="btn btn-primary">Create</button>
			</div>
		</div>
	</div>
</div>

<!-- Comment Form Dialog BS -->
<div class="modal" id="cmdlg" tabindex="-1" role="dialog" aria-labelledby="cmDLGLabel">
	<div class="modal-dialog modal-sm" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="cmDLGLabel">View Comment</h4>
			</div>
			<div class="modal-body">
				<div class="form-group">
					<label for="cmdlg_position">Position:</label>
					<input readonly type="text" id="cmdlg_position" style="background-color: #FFF;" class="form-control">
				</div>
				<div class="form-group">
					<label for="cmdlg_author">Author:</label>
					<input readonly type="text" id="cmdlg_author" style="background-color: #FFF;" class="form-control">
				</div>
				<div class="form-group">
					<label for="cmdlg_comment">Comment:</label>
					<textarea readonly id="cmdlg_comment" style="background-color: #FFF;" class="form-control" rows="3" placeholder="No comment was given."></textarea>
				</div>
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
				<button id="cmdlg_goto_btn" type="button" class="btn btn-primary">Goto Region</button>
				<button id="cmdlg_del_btn" type="button" class="btn btn-danger">Delete</button>
			</div>
		</div>
	</div>
</div>

<!-- Config Saver Dialog BS -->
<div class="modal" id="configSaveDLG" tabindex="-1" role="dialog" aria-labelledby="configSavLabel">
	<div class="modal-dialog modal-sm" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="configSavLabel">Save Current Config</h4>
			</div>
			<div class="modal-body">
				<form>
					<div id="cfg-save-error" style="display:none" class="hidden"></div>
					<div class="form-group">
						<label for="cfg_save_name">Name:</label>
						<span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="The name of the saved configuration."></span>
						<input type="text" id="cfg_save_name" class="form-control">
					</div>
					<div class="form-group">
						<label for="cfg_save_description">Description:</label>
						<span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="Short description to help identify the configuration."></span>
						<textarea id="cfg_save_description" class="form-control" rows="3" placeholder="No comment was given."></textarea>
					</div>
				</form>
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
				<button id="cfg_save_btn" type="button" class="btn btn-primary">Save</button>
			</div>
		</div>
	</div>
</div>

<!-- Config Loader Dialog BS -->
<div class="modal" id="configLoadDLG" tabindex="-1" role="dialog" aria-labelledby="configLoadLabel">
	<div class="modal-dialog" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="configLoadLabel">Load Saved Config</h4>
			</div>
			<div class="modal-body">
				<p>Manage your saved configurations here:</p>
				<table id="savedConfigTable" style="width:100%;" class="table table-striped table-bordered table-hover">
					<thead>
						<tr><th>Name</th><th>Description</th><th>Genome</th><th>Time</th><th>_TIERS</th></tr>
					</thead>
					<tbody>
					
					</tbody>
				</table>
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
				<button id="cfg_save_btn" type="button" class="btn btn-primary">Save</button>
			</div>
		</div>
	</div>
</div>

<!-- Config Manager Dialog BS -->
<div class="modal" id="configManagerDLG" style="z-index: 10001 !important;" tabindex="-1" role="dialog" aria-labelledby="cfg-manager-label">
	<div class="modal-dialog modal-sm" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="cfg-manager-label">Config</h4>
			</div>
			<div class="modal-body">
				<div class="form-group">
					<button id="cfgman-load" class="btn btn-default btn-block">Load Config</button>
				</div>
				<div class="form-group">
					<button id="cfgman-delete" class="btn btn-danger btn-block">Delete Config</button>
				</div>
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>
			</div>
		</div>
	</div>
</div>

<!-- Generic Genome Creator BS -->
<div class="modal" id="ggcDLG" style="z-index: 9099 !important;" tabindex="-1" role="dialog" aria-labelledby="ggcDLGlabel">
	<div class="modal-dialog" role="document">
		<div class="modal-content">
			<div class="modal-header">
				<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
				<h4 class="modal-title" id="ggcDLGlabel">Create Genome</h4>
			</div>
			<div class="modal-body">
				<form>
					<p>This form will allow you to create a custom genome for any species. A sequence file in 2bit form is required but the
					gene annotation file is optional. See the <a href="#" data-toggle="modal" data-target="#filemanagerDLG">file manager</a>
					to find the IDs of the required files. </p>
					<div id="ggc_form_error" style="display:none" class="hidden"></div>
					<div class="form-group">
						<label for="ggc_genome_name">Genome Name:</label>
						<span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="REQUIRED: The name of the new Genome."></span>
						<input type="text" id="ggc_genome_name" class="form-control">
					</div>
					<div class="form-group">
						<label for="ggc_seq_file">Sequence File ID:</label>
						<span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="REQUIRED: The ID of the 2bit file to use as the sequence data, as seen in the File Manager."></span>
						<input type="text" id="ggc_seq_file" class="form-control">
					</div>
					<div class="form-group">
						<label for="ggc_gene_file">Gene File ID:</label>
						<span style="margin-left: 6px;" class="glyphicon glyphicon-question-sign" aria-hidden="true" data-toggle="tooltip" data-placement="bottom" title="OPTIONAL: gene annotation file ID as seen in the file manager."></span>
						<input type="text" id="ggc_gene_file" class="form-control">
					</div>
					<input type="submit" tabindex="-1" class="hidden">
				</form>
			</div>
			<div class="modal-footer">
				<button type="button" class="btn btn-default" data-dismiss="modal">Cancel</button>
				<button id="ggc_create_btn" type="button" class="btn btn-primary">Create</button>
			</div>
		</div>
	</div>
</div>