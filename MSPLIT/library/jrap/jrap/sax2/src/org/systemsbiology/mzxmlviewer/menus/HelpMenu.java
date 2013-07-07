/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as        *
 *   published by the Free Software Foundation; either version 2 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 ***************************************************************************/

// ==================================================================
// FILE:		HelpMenu.java
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
public class HelpMenu extends JMenu
{
    private static HelpMenu _instance = new HelpMenu();
    
    private HelpMenu()
    {
        super("Help");
    }
    
    public static HelpMenu getInstance()
    {
        return _instance;
    }
}
