<html lang="en">

<head>

	<meta charset="utf-8">
	<meta http-equiv="X-UA-Compatible" content="IE=edge">
	<meta name="viewport" content="width=device-width, initial-scale=1">
	<meta name="description" content="">
	<meta name="author" content="PS">

	<title>PBrowse</title>

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
				<a id="pm_toggle" class="navbar-brand" href="/">PBrowse <span id="loggedInUser"></span></a>
			</div>
			<!-- Collect the nav links, forms, and other content for toggling -->
			<div class="collapse navbar-collapse" id="bs-example-navbar-collapse-1">
				<ul class="nav navbar-nav">
				</ul>
			</div>
			<!-- /.navbar-collapse -->
		</div>
		<!-- /.container -->
	</nav>

	<!-- Page Content -->
	<div class="container">

		<div class="row">

			<div class="col-md-2">
			</div>

			<div class="col-md-8">
				<%
				String s = (String) request.getAttribute("status");
				if ( s.equals("success") == true ) {
				%>
				<h1>Authentication Success</h1>
				<% 
				} else {
				%>		
				<h1>Authentication Failed</h1>
				<%}%>			
			
				<p><%= request.getAttribute("error") %> </p>
			</div>
			
			<div class="col-md-2">
			</div>
			
		</div>

	</div>
	<!-- /.container -->

	<div class="container">
		<hr>
		<!-- Footer -->
		<footer>
			<div class="row">
				<div class="col-lg-12">
					<p>PBrowse 2015, University of Sydney, in collaboration with VCCRI.</p>
				</div>
			</div>
		</footer>
	</div>
	<!-- /.container -->

</body>

</html>
