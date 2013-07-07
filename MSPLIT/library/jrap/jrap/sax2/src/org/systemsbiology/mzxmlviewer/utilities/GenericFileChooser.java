/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as        *
 *   published by the Free Software Foundation; either version 2 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 ***************************************************************************/

// ==================================================================
// FILE:		GenericFileChooser.java
// DATE:		May 4, 2004
// AUTHOR:		Adam Jordens (adam@jordens.org)
//
// $Id$
// Copyright 2003-2004 Adam Jordens.  All Rights Reserved.
// ==================================================================

package org.systemsbiology.mzxmlviewer.utilities;

import java.io.File;

import javax.swing.JFileChooser;
import javax.swing.filechooser.FileFilter;

/**
 * A generic file chooser that contains 1..* description/extension combinations.
 * 
 * @author <a href="mailto:adam@jordens.org">Adam Jordens</a>
 * @version $Id$
 */
public class GenericFileChooser extends JFileChooser
{
    private static GenericFileChooser _instance = new GenericFileChooser();
    
    /**
     * Default constructor.
     */
    private GenericFileChooser() 
    {
        initialize();
    }

    /**
     * @return Current instance of the GenericFileChooser.
     */
    public static GenericFileChooser getInstance()
    {
        return _instance;
    }
    
    /**
     * Add a new choosable item (description->extensions) filter.
     * 
     * @param extensions	Valid file type extensions
     * @param description	File type description
     */
    public void addFileFilter(String[] extensions, String description)
    {
        addChoosableFileFilter(new GenericFilter(extensions, description));
    }
    
    /**
     * Initialize this file chooser.
     */
    private void initialize()
    {
//        setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
    }

    /**
     * A custom FileFilter that encapsulates the accept logic with a list of
     * valid file type extensions.
     * 
     * @author <a href="mailto:adam@jordens.org">Adam Jordens</a>
     * @version $Id$
     */
    private class GenericFilter extends FileFilter
    {
        private String[] extensions;
        private String description;
        
        /**
         * Default constructor.
         * 
         * @param extensions	Valid file type extensions.
         * @param description	File type description.
         */
        public GenericFilter(String[] extensions, String description)
        {
            this.extensions = extensions;
            this.description = description;
        }
        
        /**
         * @see javax.swing.filechooser.FileFilter#accept(java.io.File)
         */
        public boolean accept(File f) 
        {
            String ext = getExtension(f);
            for (int i=0; i<extensions.length; i++)
            {
                if (extensions[i].equalsIgnoreCase(ext))
                    return true;
            }
            
            if (f.isDirectory())
                return true;
            
            return false;
        }

        /**
         * @see javax.swing.filechooser.FileFilter#getDescription()
         */
        public String getDescription() 
        {
            return description;
        }

        /**
         * @param f	File object
         * @return the extension of the File parameter.
         */
        private String getExtension(File f) 
        {
            String ext = null;
            String s = f.getName();
            int i = s.lastIndexOf('.');

            if (i > 0 &&  i < s.length() - 1) {
                ext = s.substring(i+1).toLowerCase();
            }
            return ext;
        }        
    }    
}
