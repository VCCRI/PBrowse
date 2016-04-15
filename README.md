# PBrowse: A web-based platform for real-time collaborative exploration and sharing of genomic data

Authors: Peter S. Szot, Andrian Yang, Joshua Ho

Contact: j.ho@victorchang.edu.au

Copyright (c) 2016, Victor Chang Cardiac Research Institute

This repository contains both the PBrowse web-application and the PBrowse collaborative server.

## Requirements

The server was compiled with *Java8* using *Maven* and was deployed on an *Jetty9* instance. A MySQL server was also used as the backend database software. To ensure workability, use these software packages.

```
$ apt-get install mysql-server maven
```

Java8 installation may vary depending on your OS distribution. Jetty 9 can be installed by downloading the following archive and extracting it to */opt/jetty/*:

```
$ wget http://download.eclipse.org/jetty/stable-9/dist/jetty-distribution-9.3.6.v20151106.tar.gz
```

Be sure to create the directory:

```
$ mkdir /opt/pbrowse/
```

Grant read and write access to the user group which controls the webserver. User uploaded files are stored here and much functionality will be broken without this access.The PBrowse server will create any needed subdirectories as they are needed.

The web-application was developed for the Mozilla Firefox, Chrome, and Safari browsers. IE and older browsers may have unexpected issues so use them at your own risk.

## Building

The repository should first be cloned:

```
$ git clone https://petszo@git.victorchang.edu.au/scm/holab/pbrowse.git
```

Initialise the MySQL server, as the super-user run the commands provided in the *pbrowsedb.sql* file to setup the required tables and users with either of the following commands:

```
$ mysql < pbrowsedb.sql
```
OR
```
mysql> source pbrowsedb.sql
```

All further database access is made by the application through the *pbrowse* user.

Then navigate into the *pbrowse* directory. Run the command:

```
$ mvn clean install
```

Copy the created *.war* file from *./target/* directory to your jetty installation *webapps* directory, renaming the output *war* file to *pbrowse.war* e.g.
```
$ cp ./target/*.war /opt/jetty/webapps/root.war
```

## Configure Jetty SSL

The application itself can run over http, but the underlying websocket connection must be run over SSL so you have to configure the SSL conector in Jetty. If you aren't familiar with the process, see this
[guide](http://examples.javacodegeeks.com/enterprise-java/jetty/jetty-ssl-configuration-example/) to generating a self-signed SSL certificate for testing purposes. Now add the following lines into
your */opt/jetty/etc/jetty-ssl-context.xml* file, updating the *KeyStorePath* path and password fields if you have changed them:

```
<Set name="useCipherSuitesOrder"><Property name="jetty.sslContext.useCipherSuitesOrder" default="true"/></Set>
<Set name="KeyStorePath">/opt/pbrowse/.keystore</Set>
<Set name="KeyStorePassword">pbrowse</Set>
<Set name="KeyManagerPassword">pbrowse</Set>
<Set name="TrustStorePath">/opt/pbrowse/.keystore</Set>
<Set name="TrustStorePassword">pbrowse</Set>
```

Also modify your */opt/jetty/start.ini* file and add the following two lines to the end:

```
--module=ssl
--module=https
```

## Running

Start the jetty server:

```
$ cd /opt/jetty/
$ java -jar start.jar
```

and begin browsing by visiting:
```
http://localhost:8080/
```
