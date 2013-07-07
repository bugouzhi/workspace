/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as        *
 *   published by the Free Software Foundation; either version 2 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 ***************************************************************************/

// ==================================================================
// FILE:		FileMenu.java
// DATE:		Jan 17, 2005
// AUTHOR:		Adam Jordens (adam@genologics.com)
// VERSION:     $Id$
// 
// Copyright 2005 GenoLogics Life Sciences Software, Inc.  All Rights Reserved.
// ==================================================================

package org.systemsbiology.mzxmlviewer.menus;

import javax.swing.JMenu;


/**
 * @author Adam Jordens
 */
public class FileMenu extends JMenu
{
    private static FileMenu _instance = new FileMenu();
    
    private FileMenu()
    {
        super("File");
    }
    
    public static FileMenu getInstance()
    {
        return _instance;
    }
}
