<!-- Testing interface - only written given specified parameter -->
<div id="questions">
	<div id="q0" class="test_q">
		<p>Question 0 - Me and my group</p>
		<p>
			<label>Please form yourselves into groups of 3-4 members for the following set of questions. <br><br>
			Enter your name (same name as used in solo-tests):</label>
			<input type="text" id="username" class="text ui-widget-content ui-corner-all"></input>
		</p>
	</div>
	<div id="q1" class="test_q">
		<p>Group Setup:</p>
		<p>
			<label>
			Enter a group name here which must be common to all your members: </label>
			<input type="text" id="testgroupname" class="text ui-widget-content ui-corner-all"></input>
		</p>
	</div>
	<div id="q2" class="test_q">
		<p>Are you group Leader? There can only be <em>1 per group!</em></p>
		<p>
			<label>Enter "yes" or "no":</label>
			<input type="text" id="groupleaderstatus" class="text ui-widget-content ui-corner-all"></input>
		</p>
	</div>
	<div id="q3" class="test_q">
		<p>Collaborative Communication - <span class="t_browser_name">UCSC</span></p>
		<div class="t_leaderquestion">
			<p>
			RULES:<br>
			 - You are not permitted to directly tell your followers the genes of interest. <br>
			 - Communication is allowed only via either the Genome Browser, VICTORCHANG EMAIL <br>
			 - Inform your followers of when the task has began, they will start their timers <br>
			TASK: <br>
			Have your followers visit the following genomic windows and note down all the REFSEQ genes within them: 
			</p>
			<ul class="questionlist">
				<li>chr20:1,949,623-1,974,146</li>
				<li>chr20:35,325,692-35,437,836</li>
				<li>chr20:53,177,783-53,641,205</li>
				<li>chr10:13,617,932-14,357,212</li>
				<li>chr10:117,075,807-117,815,087</li>
			</ul>
		</div>
		<div class="followerquestion">
			<p>Provided genomic locations by your group leader, note which genes are present within each coordinate window.<br>
			</p>
		</div>
			<p>Answer below:</p>
			<div class="answerField">
				<input type="text" id="genes1-UCSC" class="text ui-widget-content ui-corner-all"></input>
				<input type="text" id="genes2-UCSC" class="text ui-widget-content ui-corner-all"></input>
				<input type="text" id="genes3-UCSC" class="text ui-widget-content ui-corner-all"></input>
				<input type="text" id="genes4-UCSC" class="text ui-widget-content ui-corner-all"></input>
				<input type="text" id="genes5-UCSC" class="text ui-widget-content ui-corner-all"></input>
			</div>
			<button class="endResponse" style="width:100%;">Click here to lock in answers.</button>
	</div>
	<div id="q4" class="test_q">
		<p>Collaborative Communication - <span class="t_browser_name">IGV</span></p>
		<div class="t_leaderquestion">
			<p>
			RULES:<br>
			 - You are not permitted to directly tell your followers the genes of interest. <br>
			 - Communication is allowed only via either the Genome Browser, VICTORCHANG EMAIL <br>
			 - Inform your followers of when the task has began, they will start their timers <br>
			TASK: <br>
			Have your followers visit the following genomic windows and note down all the REFSEQ genes within them: 
			</p>
			<ul class="questionlist">
				<li>chr7:78,939,124-78,939,532</li>
				<li>chr7:23,600,869-23,601,277</li>
				<li>chr7:100,136,425-100,139,699</li>
				<li>chr7:135,367,375-135,370,648</li>
				<li>chr7:1,852,734-1,856,007</li>
			</ul>
		</div>
		<div class="followerquestion">
			<p>Provided genomic locations by your group leader, note which genes are present within each coordinate window.<br>
			</p>
		</div>
			<p>Answer below:</p>
			<div class="answerField">
				<input type="text" id="genes1-igv" class="text ui-widget-content ui-corner-all"></input>
				<input type="text" id="genes2-igv" class="text ui-widget-content ui-corner-all"></input>
				<input type="text" id="genes3-igv" class="text ui-widget-content ui-corner-all"></input>
				<input type="text" id="genes4-igv" class="text ui-widget-content ui-corner-all"></input>
				<input type="text" id="genes5-igv" class="text ui-widget-content ui-corner-all"></input>
			</div>
			<button class="endResponse" style="width:100%;">Click here to lock in answers.</button>
	</div>
	<div id="q5" class="test_q">
		<p>Collaborative Communication - <span class="t_browser_name">PBrowse</span></p>
		<p>
		Be logged into pbrowse, either make a new account or use "user" followed by 1-99, e.g. user7. The password will be the same: e.g. user7
		<img src="http://i.imgur.com/rIGDKK2.png">
		</p>
		<div class="t_leaderquestion">
			<p>
			Create a new Session using the left pannel. Instruct your followers to join it by providing them the SID.<br><br>
			RULES:<br>
			 - You are not permitted to directly tell your followers the genes of interest. <br>
			 - Communication is allowed only via either the Genome Browser, VICTORCHANG EMAIL <br>
			 - Inform your followers of when the task has began, they will start their timers <br>
			TASK: <br>
			
			Use your collaborative view to guide your followers to the targets, be sure to communicate to them when ready to proceed to new targets!<br><br>
			
			Have your followers visit the following genomic windows and note down all the REFSEQ genes within them: 
			</p>
			<ul class="questionlist">
				<li>10:31,601,746..31,644,932</li>
				<li>10:31,341,259..31,384,445</li>
				<li>10:30,973,379..31,016,565</li>
				<li>10:30,044,537..30,087,723</li>
				<li>9:29,804,614..29,847,800</li>
			</ul>
		</div>
		<div class="followerquestion">
			<p>
			Join the session provided by the leader, using the left panel, entering the given SID and clicking join<br>
			<img src="http://i.imgur.com/qX3eB5j.png"><br>
			Your leader will control the session by changing his view, you do not need to move anything.<br><br>
			Provided genomic locations by your group leader, note which genes are present within each coordinate window.<br>
			</p>
		</div>
			<p>Answer below:</p>
			<div class="answerField">
				<input type="text" id="genes1-pbrowse" class="text ui-widget-content ui-corner-all"></input>
				<input type="text" id="genes2-pbrowse" class="text ui-widget-content ui-corner-all"></input>
				<input type="text" id="genes3-pbrowse" class="text ui-widget-content ui-corner-all"></input>
				<input type="text" id="genes4-pbrowse" class="text ui-widget-content ui-corner-all"></input>
				<input type="text" id="genes5-pbrowse" class="text ui-widget-content ui-corner-all"></input>
			</div>
			<button class="endResponse" style="width:100%;">Click here to lock in answers.</button>
	</div>
	<div id="q6" class="test_q">
		<p>Collaborative File Sharing - <span class="t_browser_name">UCSC</span></p>
		<div class="t_leaderquestion">
			<p>
			RULES:<br>
			 - Share these locations in the most efficient means enabled by your genome browser.<br>
			 - Communication is allowed only via either the Genome Browser, VICTORCHANG EMAIL <br>
			 - Inform your followers of when the task has began, they will start their timers <br>
			TASK: <br>
			
			Copy this url: http://ec2-54-148-150-45.us-west-2.compute.amazonaws.com:8080/pbrowse/patient_reads_sorted.bam <br>
			Add it as a track in ucsc, instruct your followers to do the same. (UCSC does not permit uploads)
			<br><br>
			
			Communicate to your followers that they should navigate to the following GENE: MYBPC3 where everyone must work together to identify a window containing a: (start question to reveal)							
			</p>
			<ul class="questionlist">
				<li>mismatch</li>
			</ul>
		</div>
		<div class="followerquestion">
			<p>Your leader will share with you some bam data. He will also provide a task:<br>
			</p>
		</div>
		<p>Answer below:</p>
		<div class="answerField">
			<input type="text" id="mismatch-ucsc" class="text ui-widget-content ui-corner-all"></input>
		</div>
		<button class="endResponse" style="width:100%;">Click here to lock in answers.</button>
	</div>
	<div id="q7" class="test_q">
		<p>Collaborative File Sharing - <span class="t_browser_name">IGV</span></p>
		<div class="t_leaderquestion">
			<p>
			RULES:<br>
			 - Share these locations in the most efficient means enabled by your genome browser.<br>
			 - Communication is allowed only via either the Genome Browser, VICTORCHANG EMAIL <br>
			 - Inform your followers of when the task has began, they will start their timers <br>
			TASK: <br>
			
			Download these 2 files (same for next 2 questions): 
			<a href="http://ec2-54-148-150-45.us-west-2.compute.amazonaws.com:8080/pbrowse/files/12"> file1</a>, 
			<a href="http://ec2-54-148-150-45.us-west-2.compute.amazonaws.com:8080/pbrowse/files/13"> file2</a><br><br>
			For the purpose of this experiment you must pretend these files are not already hosted elsewhere and that your group members do not have copies of these files.<br>
			Share them with all the members of your group via email attachment (they are small). <br><br>
			
			Communicate to your followers that they should navigate to the following GENE: MYBPC3 where everyone must work together to identify a window containing a: (start question to reveal)							
			</p>
			<ul class="questionlist">
				<li>"G" insertion</li>
			</ul>
		</div>
		<div class="followerquestion">
			<p>Your leader will share with you some bam data. He will also provide a task:<br>
			</p>
		</div>
		<p>Answer below:</p>
		<div class="answerField">
			<input type="text" id="igv-insertion" class="text ui-widget-content ui-corner-all"></input>
		</div>
		<button class="endResponse" style="width:100%;">Click here to lock in answers.</button>
	</div>
<div id="q8" class="test_q">
	<p>Collaborative File Sharing - <span class="t_browser_name">PBrowse</span></p>
	<p>
		Be logged into pbrowse, either make a new account or use "user" followed by 1-99, e.g. user7. The password will be the same: e.g. user7
		<img src="http://i.imgur.com/rIGDKK2.png">
	</p>
	<div class="t_leaderquestion">
		<p>
		Create a new Session using the left pannel. Instruct your followers to join it by providing them the SID.<br><br>
		
		RULES:<br>
		 - Share these locations in the most efficient means enabled by your genome browser.<br>
		 - Communication is allowed only via either the Genome Browser, VICTORCHANG EMAIL <br>
		 - Inform your followers of when the task has began, they will start their timers <br>
		TASK: <br>
		
		Download these 2 files (same for next 2 questions): 
		<a href="http://ec2-54-148-150-45.us-west-2.compute.amazonaws.com:8080/pbrowse/files/12"> file1</a>, 
		<a href="http://ec2-54-148-150-45.us-west-2.compute.amazonaws.com:8080/pbrowse/files/13"> file2</a><br><br>
		For the purpose of this experiment you must pretend these files are not already hosted elsewhere and that your group members do not have copies of these files.<br>
		Share them with all the members of your group by first UPLOADING them to PBROWSE, and ADDING THEM AS A TRACK. Users in your session will share files with you automatically.<br><br>
		
		Direct your followers to the following GENE: MYBPC3, AND UNCHECK "FOLLOWER SYNC" in the LEFT PANNEL to allow followers their own exploration.<br>
		Everyone must work together to identify a window containing a single: (start question to reveal)							
		</p>
		<ul class="questionlist">
			<li>"T" insertion</li>
		</ul>
	</div>
	<div class="followerquestion">
		<p>
		Join the session provided by the leader, using the left panel, entering the given SID and clicking join<br>
		<img src="http://i.imgur.com/qX3eB5j.png"><br>
		Your leader will share with you some bam data. He will also provide a task:<br>
		Your leader will allow control of your own session, use chat and panning to locate the target and inform your group. 
		<br>
		</p>
	</div>
	<p>Answer below:</p>
	<div class="answerField">
		<input type="text" id="pbrowse-insertion" class="text ui-widget-content ui-corner-all"></input>
	</div>
	<button class="endResponse" style="width:100%;">Click here to lock in answers.</button>
</div>
<button class="startResponse hidden" style="width:100%;">Show task (Starts a timer)</button>
<button id="nextQuestion" style="width:100%;">Next Question</button>
</div>
<script>
var currentQuestion = 0;
var allResponses = {"test":"collab"};

$(function() 
{
	var tester = new Tester();
	
	$("#q0").show();
	
	var startResp = null;
	var endResp = null;
	
	//goto next question
	$("#nextQuestion").on("click",function()
	{
		var resp = {};
		var answers = $("#q"+currentQuestion).find("input");
		for (var i=0; i<answers.length; i++)
		{
			resp[answers[i].id] = answers[i].value;	
		}
		
		allResponses["q"+currentQuestion] = resp;
		
		//goto next question
		currentQuestion++;
		$(".question").hide();
		
		if ($("#q"+currentQuestion).length == 0)
		{
			//reached end
			$("#nextQuestion").hide();
			
			allResponses.timestamp = new Date().toString();
			
			var posting = $.post( 'testlogger', { result: JSON.stringify(allResponses) } );
			
			posting.done(function( data ) 
			{
				alert(data);
			});
			
			return;
		}
		$("#q"+currentQuestion).show();
		
		//hide next question if there is a timed response prompt
		var t = $("#q"+currentQuestion).find("button");
		if (t.length > 0)
		{
			$(".startResponse").addClass("hidden");
			$("#nextQuestion").hide();
		}
		$(".questionlist").hide();
		
		if (allResponses.q2 == undefined)
			return;
		
		if (allResponses.q2.groupleaderstatus == "yes")
		{
			$(".t_leaderquestion").show();
			$(".followerquestion").hide();
			$(".startResponse").removeClass("hidden");
		}
		else
		{
			$(".t_leaderquestion").hide();
			$(".startResponse").addClass("hidden");
		}
		
		tester.setGroup(allResponses.q1.testgroupname, (allResponses.q2.groupleaderstatus == "yes")?true:false)
	});
	
	$(".startResponse").on("click",function()
	{
		$("#nextQuestion").hide();
		$(".startResponse").addClass("hidden");
		$(".answerField").show();
		$(".endResponse").show();
		$(".questionlist").show();
		
		if (allResponses.q2.groupleaderstatus == "yes")
		{
			//leader initiates question start, users lock in responses
			tester.beginQuestion(currentQuestion);
		}
		
		startResp = new Date();
	});
	
	$(".endResponse").on("click",function()
	{
		$("#nextQuestion").show();
		$(".answerField").hide();
		$(".answerField").hide();
		$(".endResponse").hide();
		$(".questionlist").hide();
		endResp = new Date();
		
		if (startResp == undefined)
			startResp = new Date();
			
		allResponses["q"+currentQuestion+"_time"] = (endResp-startResp)/1000;
		
		$("#nextQuestion").click();
	});
});
</script>
