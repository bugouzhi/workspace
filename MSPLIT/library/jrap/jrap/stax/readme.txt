This is the newly designed JRAP and its implementation using StAX parser included in jdk1.6.0.
The API doesn't change vs. the old JRAP implementation.
The javadoc included for your reference.

NOTE: you have to have jre1.6 or jdk1.6 to use the software.
//July 2008
(1)add support for mzML. Now JRAP supports both mzXML and mzML formats
(2)API was slightly changed. Peak list was of double type only no matter
what underlying representation (32-bit or 64-bit precision).
