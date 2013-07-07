/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as        *
 *   published by the Free Software Foundation; either version 2 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 ***************************************************************************/

// ==================================================================
// FILE:		MenuBar.java
// DATE:		Jan 17, 2005
// AUTHOR:		Adam Jordens (adam@genologics.com)
// VERSION:     $Id$
// 
// Copyright 2005 GenoLogics Life Sciences Software, Inc.  All Rights Reserved.
// ==================================================================

package org.systemsbiology.mzxmlviewer.menus;

import javax.swing.JMenuBar;


/**
 * @author Adam Jordens
 */
public class RootMenuBar extends JMenuBar
{
    private static RootMenuBar _instance = new RootMenuBar();
    
    private RootMenuBar()
    {
        add(FileMenu.getInstance());
        add(HelpMenu.getInstance());
    }
    
    public static RootMenuBar getInstance()
    {
        return _instance;
    }
}
