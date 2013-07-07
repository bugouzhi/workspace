/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as        *
 *   published by the Free Software Foundation; either version 2 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 ***************************************************************************/

/*******************************************************************************
 * --------------------------------------------------------------------------- *
 * File: * @(#) SAX2ScanHandler.java * Author: * Robert M. Hubley
 * rhubley@systemsbiology.org
 * ****************************************************************************** * * *
 * This software is provided ``AS IS'' and any express or implied * *
 * warranties, including, but not limited to, the implied warranties of * *
 * merchantability and fitness for a particular purpose, are disclaimed. * * In
 * no event shall the authors or the Institute for Systems Biology * * liable
 * for any direct, indirect, incidental, special, exemplary, or * *
 * consequential damages (including, but not limited to, procurement of * *
 * substitute goods or services; loss of use, data, or profits; or * * business
 * interruption) however caused and on any theory of liability, * * whether in
 * contract, strict liability, or tort (including negligence * * or otherwise)
 * arising in any way out of the use of this software, even * * if advised of
 * the possibility of such damage. * * *
 * ******************************************************************************
 * 
 * ChangeLog
 * 
 * $Log: SAX2ScanHandler.java,v $ Revision 1.1.1.1 2003/04/09 00:02:54 ppatrick
 * Initial import.
 * 
 * 10-05-2004: fixed bug in for loop, M. Vogelzang
 * 
 * 1.2 added logging & handling of precursorMZ 
 * M. Vogelzang  
 * (logging commented out to prevent need of log4j library)
 * 
 ******************************************************************************/
package org.systemsbiology.jrap;
//import org.apache.log4j.Logger;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;
import org.xml.sax.helpers.DefaultHandler;
import java.util.zip.*;

public final class SAX2ScanHandler extends DefaultHandler
{
//	private static final Logger logger =
//		Logger.getLogger(SAX2ScanHandler.class);
	/** A new Scan Object */
	protected Scan tmpScan;

	/** A string to hold the Base64 peak data */
	protected StringBuffer peakData = new StringBuffer();

	/** Flag to indicate we are reading a peak tag */
	protected boolean inPeak = false;

	/** Buffer to hold characters while getting precursorMZ value */
	protected StringBuffer precursorBuffer;

	/** Flag to indicate if we are reading the precursor MZ value */
	protected boolean inPrecursorMZ = false;

	//
	// Getters
	//

	public Scan getScan()
	{
		return (tmpScan);
	}

	private int getIntAttribute(Attributes attrs, String name)
	{
		int result;

		if (attrs.getValue(name) == null) // attribute not present
			return -1;

		try
		{
			result = Integer.parseInt(attrs.getValue(name));
		} catch (NumberFormatException e)
		{
//			logger.error("Numberformatexception!", e);
			result = -1;
		}
		return (result);
	}

	private long getLongAttribute(Attributes attrs, String name)
	{
		long result;

		if (attrs.getValue(name) == null) // attribute not present
			return -1;

		try
		{
			result = Long.parseLong(attrs.getValue(name));
		} catch (NumberFormatException e)
		{
//			logger.error("Numberformatexception!", e);
			result = -1;
		}
		return (result);
	}

	private float getFloatAttribute(Attributes attrs, String name)
	{
		float result;

		if (attrs.getValue(name) == null) // attribute not present
			return -1;

		try
		{
			result = Float.parseFloat(attrs.getValue(name));
		} catch (NumberFormatException e)
		{
//			logger.error("Numberformatexception!", e);
			result = -1;
		} catch (NullPointerException e1)
		{
//			logger.error("Nullpointerexception!", e1);
			result = -1;
		}
		return (result);
	}

    private void calFloatList(byte[] tmpArr)
    {	
	
	int floatBytes = tmpScan.getPrecision() / 8;
	float[][] tmpMassIntensityList = new float[2][tmpArr.length / floatBytes / 2];
	int peakIndex = 0;
	int fieldIndex = 0;

	if (floatBytes <= 0)
	    System.err.println("FLOATBYTES <= 0!!!");

	for (int i = 0; i < tmpArr.length - floatBytes; i += floatBytes)
	    {
		int intBits = 0;
		intBits |= (((int) tmpArr[i]) & 0xff);
		intBits <<= 8;
		intBits |= (((int) tmpArr[i + 1]) & 0xff);
		intBits <<= 8;
		intBits |= (((int) tmpArr[i + 2]) & 0xff);
		intBits <<= 8;
		intBits |= (((int) tmpArr[i + 3]) & 0xff);
		// Must be in IEEE 754 encoding!
		tmpMassIntensityList[fieldIndex++][peakIndex] =
		    Float.intBitsToFloat(intBits);
		if (fieldIndex == 2)
		    {
			fieldIndex = 0;
			peakIndex++;
		    }
	    }
	inPeak = false;
	peakData.delete(0, peakData.capacity());
	tmpScan.setMassIntensityList(tmpMassIntensityList);
	tmpMassIntensityList = null;
	tmpArr = null;
    }
			
		
    
	//
	// ContentHandler Methods
	//

	/** Start document. */
	public void startDocument() throws SAXException
	{
		// Nothing to do
	} // startDocument()

	/** Start element. */
	public void startElement(
		String uri,
		String local,
		String raw,
		Attributes attrs)
		throws SAXException
	{
		if (raw.equals("scan"))
		{
			tmpScan = new Scan();
			tmpScan.setNum(getIntAttribute(attrs, "num"));
			tmpScan.setMsLevel(getIntAttribute(attrs, "msLevel"));
			tmpScan.setPeaksCount(getIntAttribute(attrs, "peaksCount"));
			tmpScan.setPolarity(attrs.getValue("polarity"));
			tmpScan.setScanType(attrs.getValue("scanType"));
			tmpScan.setCentroided(getIntAttribute(attrs, "centroided"));
			tmpScan.setDeisotoped(getIntAttribute(attrs, "deisotoped"));
			tmpScan.setChargeDeconvoluted(
				getIntAttribute(attrs, "chargeDeconvoluted"));
			tmpScan.setRetentionTime(attrs.getValue("retentionTime"));
			tmpScan.setStartMz(getFloatAttribute(attrs, "startMz"));
			tmpScan.setEndMz(getFloatAttribute(attrs, "endMz"));
			tmpScan.setLowMz(getFloatAttribute(attrs, "lowMz"));
			tmpScan.setHighMz(getFloatAttribute(attrs, "highMz"));
			tmpScan.setBasePeakMz(getFloatAttribute(attrs, "basePeakMz"));
			tmpScan.setBasePeakIntensity(
				getFloatAttribute(attrs, "basePeakIntensity"));
			tmpScan.setTotIonCurrent(getFloatAttribute(attrs, "totIonCurrent"));
			tmpScan.setFilterLine(attrs.getValue("filterLine"));
		} else if (raw.equals("peaks"))
		{
			tmpScan.setPrecision(getIntAttribute(attrs, "precision"));
			tmpScan.setByteOrder(attrs.getValue("byteOrder"));
		      
			if(attrs.getValue("contentType") == null)
			    tmpScan.setContentType("none");
			else
			    tmpScan.setContentType(attrs.getValue("contentType"));

			if(attrs.getValue("compressionType") == null)
			    tmpScan.setCompressionType("none");
			else
			    tmpScan.setCompressionType(attrs.getValue("compressionType"));

			tmpScan.setCompressedLen(getIntAttribute(attrs, "compressedLen"));
			inPeak = true;
		} else if (raw.equals("precursorMz"))
		{
			tmpScan.setPrecursorScanNum(
				getIntAttribute(attrs, "precursorScanNum"));
			tmpScan.setPrecursorCharge(
				getIntAttribute(attrs, "precursorCharge"));
			tmpScan.setCollisionEnergy(
				getFloatAttribute(attrs, "collisionEnergy"));
			tmpScan.setIonisationEnergy(
				getFloatAttribute(attrs, "ionisationEnergy"));
			
			tmpScan.setPrecursorIntensity(getFloatAttribute
						      (attrs, "precursorIntensity"));

			precursorBuffer = new StringBuffer();
			inPrecursorMZ = true;
		}
	} // startElement(String,String,StringAttributes)

	public void endElement(String uri, String local, String raw)
	    throws SAXException
	{
		if (raw.equals("peaks"))
		{
		    byte[] tmpArr = Base64.decode(peakData.toString());

		    byte[] result = null;
		    boolean isCompressed = false;

		    if((tmpScan.getCompressionType()).equals("zlib"))
			isCompressed=true;
		    else
			isCompressed = false;
		   
		    if(isCompressed)
			{
			    try{
				Inflater decompresser = new Inflater();
				decompresser.setInput(tmpArr, 0, tmpScan.getCompressedLen());
				int unCompLen = (tmpScan.getPeaksCount())*(tmpScan.getPrecision()/4)+1;
				result = new byte[unCompLen];
				decompresser.inflate(result);
				decompresser.end();

			    }
			    catch(DataFormatException e)
				{
				    e.printStackTrace();
				}
			}
		    else 
			result=tmpArr;
			 
		    
		    calFloatList(result);

		    
		    throw (new SAXException("ScanEndFoundException"));
				
		} else if (raw.equals("precursorMz"))
		{
			tmpScan.setPrecursorMz(
				Float.parseFloat(precursorBuffer.toString()));
			precursorBuffer = null; // make available for garbage collection

			inPrecursorMZ = false;
		}
	} // endElement()

	/** Characters. */
	public void characters(char ch[], int start, int length)
		throws SAXException
	{
		if (inPeak)
		{
			peakData.append(ch, start, length);
		} else if (inPrecursorMZ)
		{
			precursorBuffer.append(ch, start, length);
		}
	} // characters(char[],int,int);

	/** Ignorable whitespace. */
	public void ignorableWhitespace(char ch[], int start, int length)
		throws SAXException
	{
		// Do nothing
	} // ignorableWhitespace(char[],int,int);

	/** Processing instruction. */
	public void processingInstruction(String target, String data)
		throws SAXException
	{
		// Do nothing
	} // processingInstruction(String,String)

	//
	// ErrorHandler methods
	//

	/** Warning. */
	public void warning(SAXParseException ex) throws SAXException
	{
		// Do nothing
		//printError("Warning", ex);
	} // warning(SAXParseException)

	/** Error. */
	public void error(SAXParseException ex) throws SAXException
	{
		// Do nothing
		//printError("Error", ex);
	} // error(SAXParseException)

	/** Fatal error. */
	public void fatalError(SAXParseException ex) throws SAXException
	{
		// Do nothing
		//printError("Fatal Error", ex);
	} // fatalError(SAXParseException)

	//
	// Protected methods
	//

	/** Prints the error message. */
	protected void printError(String type, SAXParseException ex)
	{
		System.err.print("[");
		System.err.print(type);
		System.err.print("] ");
		if (ex == null)
		{
			System.out.println("!!!");
		}
		String systemId = ex.getSystemId();
		if (systemId != null)
		{
			int index = systemId.lastIndexOf('/');
			if (index != -1)
				systemId = systemId.substring(index + 1);
			System.err.print(systemId);
		}
		System.err.print(':');
		System.err.print(ex.getLineNumber());
		System.err.print(':');
		System.err.print(ex.getColumnNumber());
		System.err.print(": ");
		System.err.print(ex.getMessage());
		System.err.println();
		System.err.flush();
	} // printError(String,SAXParseException)

}
