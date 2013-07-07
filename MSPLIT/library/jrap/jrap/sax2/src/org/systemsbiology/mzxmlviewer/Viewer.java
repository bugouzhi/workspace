/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as        *
 *   published by the Free Software Foundation; either version 2 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 ***************************************************************************/


// ==================================================================
// FILE:		Viewer.java
// DATE:		Jan 17, 2005
// AUTHOR:		Adam Jordens (adam@genologics.com)
// VERSION:     $Id$
// ==================================================================

package org.systemsbiology.mzxmlviewer;

import java.io.File;

import javax.swing.*;

import org.systemsbiology.mzxmlviewer.utilities.GenericFileChooser;

import com.jgoodies.looks.plastic.*;
import com.jgoodies.looks.plastic.theme.SkyBluerTahoma;



/**
 * @author Adam Jordens
 */
public class Viewer
{
    public static void main(String[] args)
    {
        SwingUtilities.invokeLater(
            new Runnable()
            {                
                /**
                 * @see java.lang.Runnable#run()
                 */
                public void run()
                {
                    setupLookAndFeel();                   
                    
                    GenericFileChooser fileChooser = GenericFileChooser.getInstance();
                    fileChooser.addFileFilter(new String[] {"xml", "mzXML"}, "mzXML File");
                    
                    int status = fileChooser.showOpenDialog(null);
                    
                    if (status == JFileChooser.APPROVE_OPTION)
                    {
                        File selectedFile = fileChooser.getSelectedFile();
                        new MzXMLViewer(selectedFile).show();
                    }
                }
            });       
    }
    
    /**
     * Configure the application's standard look & feel.
     */
    private static void setupLookAndFeel()
    {
        PlasticLookAndFeel.setMyCurrentTheme(new SkyBluerTahoma());
        try {
           UIManager.setLookAndFeel(new Plastic3DLookAndFeel());
        } catch (Exception e) {}   
    }
}
