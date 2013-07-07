/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as        *
 *   published by the Free Software Foundation; either version 2 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 ***************************************************************************/

// ==================================================================
// FILE:		LoadDataDialog.java
// DATE:		Jan 17, 2005 
// AUTHOR:		Adam Jordens (adam@genologics.com)
// VERSION:     $Id$
// ==================================================================
package org.systemsbiology.mzxmlviewer.utilities;

import java.awt.*;

import javax.swing.*;

import org.systemsbiology.mzxmlviewer.MzXMLViewer;

/**
 * @author Adam Jordens
 */
public class LoadDataDialog extends JDialog
{    
    public LoadDataDialog(MzXMLViewer parentFrame)
    {
        super(parentFrame, "Loading data...", true);

        setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
        setSize(500, 300);
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        setLocation((int) (screenSize.getWidth() / 2 - 250), (int) (screenSize
                .getHeight() / 2 - 150));

        Container splashContainer = new JPanel();
        splashContainer.setBackground(Color.ORANGE);
        getContentPane().setBackground(Color.ORANGE);
        getContentPane().setLayout(new GridBagLayout());
        getContentPane().add(splashContainer);
        splashContainer.setLayout(new BoxLayout(splashContainer,
                BoxLayout.Y_AXIS));
        JLabel title = new JLabel("mzXML viewer version " + MzXMLViewer.VERSION_STRING);
        title.setFont(title.getFont().deriveFont(24.0f));
        title.setAlignmentX(0.5f);
        JLabel copyR = new JLabel(
                "(C) 2004 by the Institute of Systems Biology (www.systemsbiology.org)");
        copyR.setAlignmentX(0.5f);
        JLabel pleaseWait = new JLabel(
                "Please wait while the file is being loaded...");
        pleaseWait.setAlignmentX(0.5f);
        pleaseWait.setFont(pleaseWait.getFont().deriveFont(10.0f));
        splashContainer.add(title);
        splashContainer.add(copyR);
        splashContainer.add(pleaseWait);   
    }
}
