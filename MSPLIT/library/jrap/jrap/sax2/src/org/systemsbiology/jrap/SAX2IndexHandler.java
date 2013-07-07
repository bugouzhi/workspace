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
 * File: * @(#) SAX2IndexHandler.java * Author: * Robert M. Hubley
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
 * $Log: SAX2IndexHandler.java,v $ Revision 1.1.1.1 2003/04/09 00:02:54
 * ppatrick Initial import.
 * 1.2 added buffering of characters Mathijs
 * 
 *  
 ******************************************************************************/
package org.systemsbiology.jrap;

import java.util.ArrayList;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;
import org.xml.sax.helpers.DefaultHandler;

public final class SAX2IndexHandler extends DefaultHandler
{

    /** Buffer for strings that get split over multiple characters() calls */
    private StringBuffer buffer;

    /** Flag to indicate if there is an index in this file */
    protected boolean foundIdxOffset = false;

    /** If an index exists this is the byte position where it occurs */
    protected long idxOffset = -1;

    /** True if we are sitting on an offset line */
    protected boolean foundScanOffset = false;

    /** Data structure to hold index table */
    ArrayList indexes = new ArrayList();

    //
    // Getters
    //

    /** Status of index offset */
    public boolean foundOffset()
    {
        return (foundIdxOffset);
    }

    /** Get byte position of the index in the file */
    public long getIndexPosition()
    {
        return (idxOffset);
    }

    /** Get scan offset */
    public long getScanOffset(int scanNumber)
    {
        if (scanNumber > 0)
        {
            return (((Long) indexes.get(scanNumber - 1)).longValue());
        } else
        {
            return (-1);
        }
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
    public void startElement(String uri, String local, String raw,
            Attributes attrs) throws SAXException
    {
        if (raw.equals("indexOffset"))
        {
            foundIdxOffset = true;
            buffer = new StringBuffer();
        } else if (raw.equals("offset"))
        {
            // TODO: Will there be anything but a scan index????
            if (attrs != null)
            {
                if (attrs.getLength() == 1)
                {
                    foundScanOffset = true;
                    buffer = new StringBuffer();
                    // For now we are going to assume that
                    // the id="" attribute is consecutive from
                    // 1-n.
                    //
                }
            }
        }
    } // startElement(String,String,StringAttributes)

    public void endElement(String uri, String local, String raw)
            throws SAXException
    {
        // System.out.println(" endof " + raw);
        if (raw.equals("offset"))
        {
            try
            {
                Long l = new Long(buffer.toString());
                indexes.add(l);
            } catch (NumberFormatException e)
            {
                System.err.println("Error: File contains an invalid offset!: "
                        + buffer.toString());
            }
            foundScanOffset = false;
        } else if (raw.equals("indexOffset"))
        {
            System.out.println("end of index, parsing " + buffer.toString());
            try
            {
                idxOffset = Long.parseLong(buffer.toString());
            } catch (NumberFormatException e)
            {
                System.err
                        .println("Error: File contains an invalid index offset!: "
                                + buffer.toString());
                idxOffset = -1;
            }
            throw new SAXException("IdxReadException");
        }
        else if(raw.equals("index"))
            throw new SAXException("IdxReadException");
    } // endElement()

    /** Characters. */
    public void characters(char ch[], int start, int length)
            throws SAXException
    {
        if (foundIdxOffset || foundScanOffset)
        {
            buffer.append(new String(ch, start, length));
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
        System.out.println("instr? " + target + " --- " + data);
        // Do nothing
    } // processingInstruction(String,String)

    //
    // ErrorHandler methods
    //

    /** Warning. */
    public void warning(SAXParseException ex) throws SAXException
    {
        // Do nothing
        // printError("Warning", ex);
    } // warning(SAXParseException)

    /** Error. */
    public void error(SAXParseException ex) throws SAXException
    {
        // Do nothing
        // printError("Error", ex);
    } // error(SAXParseException)

    /** Fatal error. */
    public void fatalError(SAXParseException ex) throws SAXException
    {
        // Do nothing
        // printError("Fatal Error", ex);
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
            if (index != -1) systemId = systemId.substring(index + 1);
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