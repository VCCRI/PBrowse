# PBrowse: The Collaborative Genome Browser

This repository contains both the PBrowse web-application and the PBrowse collaborative server.

## Requirements

The server was compiled with *Java8* using *Maven* and was deployed on an *Jetty9* instance. A MySQL server was also used as the backend database software. To ensure workability, use these software packages. Java8 installation may vary depending on your OS distribution.

Installing Java8

```
sudo yum install -y java-1.8.0-openjdk.x86_64 java-1.8.0-openjdk-devel.x86_64
sudo update-alternatives --config java #You will need to then enter a number corresponding to jdk8
```

Installing mysql

```
sudo yum install -y mysql-server
sudo chkconfig mysqld on
sudo service mysqld start
mysqladmin -u root password [your_new_pwd]
```

Installing maven [source](https://gist.github.com/sebsto/19b99f1fa1f32cae5d00)

```
sudo wget http://repos.fedorapeople.org/repos/dchen/apache-maven/epel-apache-maven.repo -O /etc/yum.repos.d/epel-apache-maven.repo
sudo sed -i s/\$releasever/6/g /etc/yum.repos.d/epel-apache-maven.repo
sudo yum install -y apache-maven
```

Jetty 9 can be installed by downloading the following archive and extracting it to */opt/jetty/*:

```
wget http://download.eclipse.org/jetty/stable-9/dist/jetty-distribution-9.3.6.v20151106.tar.gz
tar -xzf jetty-distribution-9.3.6.v20151106.tar.gz
sudo mkdir /opt/jetty
sudo mv -r jetty-distribution-9.3.6.v20151106/* /opt/jetty/
```

Be sure to create the directory:

```
sudo mkdir /opt/pbrowse/
```

Grant read and write access to the user group which controls the webserver. User uploaded files are stored here and much functionality will be broken without this access.The PBrowse server will create any needed subdirectories as they are needed.

Create a new user for running webserver

```
sudo adduser pbrowse
sudo chown -R pbrowse:root /opt/pbrowse
sudo chown -R pbrowse:root /opt/jetty
```

The web-application was developed for the Mozilla Firefox, Chrome, and Safari browsers. IE and older browsers may have unexpected issues so use them at your own risk.

## Building

The repository should first be cloned:

```
sudo yum install -y git
git clone https://[username]@git.victorchang.edu.au/scm/holab/pbrowse.git
```

Initialise the MySQL server, as the super-user run the commands provided in the *pbrowsedb.sql* file to setup the required tables and users with either of the following commands:

```
mysql -u root -p < pbrowsedb.sql
```
OR
```
mysql> source pbrowsedb.sql
```

All further database access is made by the application through the *pbrowse* user.

Then navigate into the *pbrowse* directory. Run the command:

```
mvn clean install
```

Copy the created *.war* file from *./target/* directory to your jetty installation *webapps* directory, renaming the output *war* file to *pbrowse.war* e.g.
```
cp ./target/*.war /opt/jetty/webapps/root.war
```

## Configure Jetty SSL

The application itself can run over http, but the underlying websocket connection must be run over SSL so you have to configure the SSL connector in Jetty. If you aren't familiar with the process, see this
[guide](https://wiki.eclipse.org/Jetty/Howto/Configure_SSL#Generating_Key_Pairs_and_Certificates) to generating a self-signed SSL certificate for testing purposes.

```
keytool -genkey -keyalg RSA -alias pbrowse -keystore pbrowse.jks -validity 365 -keysize 2048
keytool -importkeystore -srckeystore pbrowse.jks -destkeystore keystore
sudo cp keystore /opt/jetty/etc
```

For production, you will need to purchase SSL from a CA (such as GeoTrust, Comodo etc). Once you have purchased, you will normally be provided with 3 files - a .key file containing your private key, a .crt file containing the certificate for your domain and another .crt file containing intermediate CA certificate (if you buy the certificate from a reseller like goDaddy and RapidSSL). The instruction below is tailored to RapidSSL specifically (Instruction taken from answer by s_t_e_v_e in [here](http://stackoverflow.com/questions/4008837/configure-ssl-on-jetty)).
```
cat your.domain.crt intermediateCA.crt  > certchain.crt # combine the two certificate to one file
openssl pkcs12 -export -inkey your.privatekey.key -in certchain.crt -out jetty.pkcs12 # create a pkcs12 file (you will be asked to set up a password for pkcs12 file)
keytool -importkeystore -srckeystore jetty.pkcs12 -srcstoretype PKCS12 -destkeystore keystore (you will need to set up another password for the keystore (which will be the password used for the settings below)
```

Check to make sure that your keystore contains your certificate chains
```
keytool -list -v -keystore keystore
```

Expected output
```
Keystore type: JKS
Keystore provider: SUN

Your keystore contains 1 entry

Alias name: 1
Creation date: 17/02/2016
Entry type: PrivateKeyEntry
Certificate chain length: 2
Certificate[1]:
...
Certificate[2]:
...
```

Now add the following lines into your */opt/jetty/etc/jetty-ssl-context.xml* file, updating the *KeyStorePath* path (currently set up */opt/jetty/etc*) and password fields if you have changed them:

```
<Set name="useCipherSuitesOrder"><Property name="jetty.sslContext.useCipherSuitesOrder" default="true"/></Set>
<Set name="KeyStorePath">/opt/jetty/etc/keystore</Set>
<Set name="KeyStorePassword">pbrowse</Set>
<Set name="KeyManagerPassword">pbrowse</Set>
<Set name="TrustStorePath">/opt/jetty/etc/keystore</Set>
<Set name="TrustStorePassword">pbrowse</Set>
```

Also modify your */opt/jetty/start.ini* file and add the following two lines to the end:

```
--module=ssl
--module=https
```

## Configure email bot

You will need to create *email_bot_config.properties* under */opt/pbrowse*. The file needs to contain

```
email=[username]@gmail.com
password=password
```

Currently only gmail account is supported as the smtp settings is hardcoded to use gmail's. This can be modified on the SendMailTLS class.

## Running

Start the jetty server:

```
sudo -su pbrowse
cd /opt/jetty/
java -jar start.jar &
```

Note: & is required to run jetty in the background

and begin browsing by visiting:
```
http://localhost:8080/
```


## Enabling access to standard ports

To access pbrowse, you currently need to specify the 8080 port for http and 8443 port for https. In order to allow access using the standard http and https port (80 and 443 respectively), you need to set up iptables redirect. (Information taken from [jetty doc](http://www.eclipse.org/jetty/documentation/current/setting-port80-access.html) and [this instruction](https://gist.github.com/kentbrew/776580).)

Check that no configuration has been previously made using 
```
sudo iptables -t nat -L
```

You should see the output below:
```
Chain PREROUTING (policy ACCEPT)
target     prot opt source               destination         

Chain INPUT (policy ACCEPT)
target     prot opt source               destination         

Chain OUTPUT (policy ACCEPT)
target     prot opt source               destination         

Chain POSTROUTING (policy ACCEPT)
target     prot opt source               destination
```

Now, set up the rules for http and https
```
sudo iptables -t nat -I PREROUTING -p tcp --dport 80 -j REDIRECT --to-port 8080
sudo iptables -t nat -I PREROUTING -p tcp --dport 443 -j REDIRECT --to-port 8443
```

If you run the command to check the configuration again (`sudo iptables -t nat -L`), you should now see the new rules
```
Chain PREROUTING (policy ACCEPT)
target     prot opt source               destination         
REDIRECT   tcp  --  anywhere             anywhere             tcp dpt:https redir ports 8443
REDIRECT   tcp  --  anywhere             anywhere             tcp dpt:http redir ports 8080

Chain INPUT (policy ACCEPT)
target     prot opt source               destination         

Chain OUTPUT (policy ACCEPT)
target     prot opt source               destination         

Chain POSTROUTING (policy ACCEPT)
target     prot opt source               destination
```