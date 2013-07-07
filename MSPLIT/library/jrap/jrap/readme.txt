the JRAP package was implemented using SAX2 parser from apache and this was included in sax2 directory.
Recently the JRAP was re-designed and re-implemented using StAX parser from jdk1.6 and this was included
in the stax directory.

In order to hide the implementation from users, the API didn't change. You don't need to change your code
if you use the new JRAP based on StAX.

NOTE: you have to have jre1.6 or jdk1.6 to use the StAX version JRAP.
