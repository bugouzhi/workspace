/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as        *
 *   published by the Free Software Foundation; either version 2 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 ***************************************************************************/

/*******************************************************************************
 * --------------------------------------------------------------------------- 
 * File: MzTableModel.java  
 * Author: Mathijs Vogelzang
 * m_v@dds.nl
 * ****************************************************************************** 
 * This software is provided ``AS IS'' and any express or implied
 * warranties, including, but not limited to, the implied warranties of 
 * merchantability and fitness for a particular purpose, are disclaimed.
 * In no event shall the authors or the Institute for Systems Biology liable
 * for any direct, indirect, incidental, special, exemplary, or
 * consequential damages (including, but not limited to, procurement of
 * substitute goods or services; loss of use, data, or profits; or business
 * interruption) however caused and on any theory of liability, whether in
 * contract, strict liability, or tort (including negligence or otherwise)
 * arising in any way out of the use of this software, even if advised of
 * the possibility of such damage.
 * *****************************************************************************
 */

package org.systemsbiology.mzxmlviewer.tables;

import javax.swing.table.AbstractTableModel;

import org.systemsbiology.jrap.ScanHeader;

/**
 * MzTableModel is a concrete implementation of AbstractTableModel
 * to make a table showing scans.
 * 
 * @author M. Vogelzang
 */
public class MzTableModel extends AbstractTableModel
{
	private ScanHeader[] scanHeaders;

	public int getColumnCount()
	{
		return 3;
	}

	public int getRowCount()
	{
		return scanHeaders.length;
	}

	public Object getValueAt(int rowIndex, int columnIndex)
	{
		switch (columnIndex)
		{
			case 0 :
				return new Integer(scanHeaders[rowIndex].getNum());
			case 1 :
				return new Integer(scanHeaders[rowIndex].getMsLevel());
			case 2 :
				if (scanHeaders[rowIndex].getMsLevel() == 1)
					return "---";
				else
					return new Double(scanHeaders[rowIndex].getPrecursorMz());
		}
		return null;
	}

	public MzTableModel(ScanHeader[] scanHeaders)
	{
		this.scanHeaders = scanHeaders;
	}

	public String getColumnName(int column)
	{
		switch (column)
		{
			case 0 :
				return "Scan number";
			case 1 :
				return "MS level";
			case 2 :
				return "Precursor mass";
		}
		return "???";
	}

}
