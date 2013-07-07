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
 * File: MzXMLViewer.java  
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
 * 
 * Changelog:
 * 10-05-2004: M. Vogelzang, adapted everything to xerces 2_6_2 and jfreechart
 * 0.9.21, added VERSION_STRING
 */
package org.systemsbiology.mzxmlviewer;

import java.awt.*;
import java.awt.event.*;
import java.io.*;

import javax.swing.*;

import org.systemsbiology.jrap.*;
import org.systemsbiology.mzxmlviewer.dialogs.InfoDialog;
import org.systemsbiology.mzxmlviewer.menus.*;
import org.systemsbiology.mzxmlviewer.tables.MzTableModel;
import org.systemsbiology.mzxmlviewer.utilities.*;

/**
 * MzXMLViewer is a simple viewer for mzXML files.
 * 
 * @author M. Vogelzang
 */
public class MzXMLViewer 
    extends JFrame
	implements ActionListener, SpectrumComponent.ClickListener
{
	public static final String VERSION_STRING = "1.4";

	private MzTableModel tableModel;
	private SpectrumComponent spectrum;
	private MSXMLParser parser;
	private ScanHeader[] scanHeaders;
	private CardLayout lowerLayout;
	private JComboBox selector;
	private JPanel cardPanel;
	private JTable table;

	private JMenuItem openMenuItem, showInfoMenuItem, exitMenuItem,
			aboutMenuItem;

	private final File mzXMLFile;

    /**
     * Construct the base MzXMLViewer Frame.
     * 
     * @param filename  The file that should be parsed and visualized.
     */
    public MzXMLViewer(final File mzXMLFile)
    {
        super("mzXML spectrum viewer version " + VERSION_STRING + " | "
                + mzXMLFile.getName());
                
        this.mzXMLFile = mzXMLFile;
        
        refreshData();
        
        setExtendedState(MAXIMIZED_BOTH);
        setDefaultCloseOperation(EXIT_ON_CLOSE);
        pack();        
    }    
        
    /**
     * Display the currently selected scan.
     * @param scanNumber    The scan to display.
     */
	public void selectScan(int scanNumber)
	{
		table.setRowSelectionInterval(scanNumber - 1, scanNumber - 1);
		table.scrollRectToVisible(new Rectangle(0, (scanNumber - 1)
				* (table.getRowHeight()), 20, 20));
		spectrum.setScan(parser.rap(scanNumber));
	}

    /**
     * Refresh/re-parse the data from the selected file.
     */
    private void refreshData()
    {
        final LoadDataDialog loadWorker = new LoadDataDialog(this);
        
        new SwingWorker()
        {            
            /**
             * @see org.systemsbiology.mzxmlviewer.utilities.SwingWorker#construct()
             */
            public Object construct()
            {
                try
                {
                    parser = new MSXMLParser(mzXMLFile.getCanonicalPath());
                }
                catch(IOException e)
                {
                    throw new IllegalStateException("Unable to parse mzXML File : " + mzXMLFile);
                }
                
                scanHeaders = new ScanHeader[parser.getScanCount()];
                for (int i=1; i<=scanHeaders.length; i++)
                {
                    scanHeaders[(i - 1)] = parser.rapHeader(i);
                }
                              
                return null;
            }
                        
            /**
             * @see org.systemsbiology.mzxmlviewer.utilities.SwingWorker#finished()
             */
            public void finished()
            {
                loadWorker.dispose();
            }
        }.start();
        
        loadWorker.show();
        
        // Remove all GUI components and re-add them.
        getContentPane().removeAll();
        setupComponents();
    }
    
    /**
     * Setup this MzXMLViewer's GUI Components.  
     */
	private void setupComponents()
	{
		Container c = getContentPane();
		JPanel c2 = new JPanel(new BorderLayout(11, 0));
		tableModel = new MzTableModel(scanHeaders);

		table = new JTable(tableModel);
		table.setPreferredScrollableViewportSize(new Dimension(400, 200));
		c2.add(new JScrollPane(table), BorderLayout.WEST);

		table.addMouseListener(new MouseAdapter()
		{
			public void mouseClicked(MouseEvent me)
			{
				if (table.getSelectedRow() != -1)
				{
					selectScan(table.getSelectedRow() + 1);
				}
			}
		});

		spectrum = new SpectrumComponent();
		c2.add(spectrum);

		JPanel lowerPanel = new JPanel();
		buildLowerPanel(lowerPanel);

		c2.setBorder(BorderFactory.createEmptyBorder(11, 11, 11, 11));
		c.add(c2, BorderLayout.CENTER);
		lowerPanel.setBorder(BorderFactory.createEmptyBorder(0, 11, 11, 11));
		c.add(lowerPanel, BorderLayout.SOUTH);

		RootMenuBar rootMenuBar = RootMenuBar.getInstance();
        FileMenu file = FileMenu.getInstance();
        HelpMenu help = HelpMenu.getInstance();

        // TODO:  Refactor these menu items into their respective menu's.
        
		openMenuItem = new JMenuItem("Open...");
		openMenuItem.addActionListener(this);
		showInfoMenuItem = new JMenuItem("Show header info");
		showInfoMenuItem.addActionListener(this);
		exitMenuItem = new JMenuItem("Exit");
		exitMenuItem.addActionListener(this);
		aboutMenuItem = new JMenuItem("About");
		aboutMenuItem.addActionListener(this);
		file.add(openMenuItem);
		file.add(showInfoMenuItem);
		file.addSeparator();
		file.add(exitMenuItem);

		help.add(aboutMenuItem);
		setJMenuBar(rootMenuBar);

		selectScan(1);
	}

	private void buildLowerPanel(JPanel panel)
	{
		panel.setLayout(new BorderLayout());

		lowerLayout = new CardLayout();
		cardPanel = new JPanel(lowerLayout);

		cardPanel.add(SpectrumComponent.getTICComponent(scanHeaders, this),
				"TIC");
		cardPanel.add(SpectrumComponent.getBPIComponent(scanHeaders, this),
				"BPI");

		JPanel upperPanel = new JPanel();
		selector = new JComboBox(new String[]
		{"Total ion current", "Base peak intensity"});
		selector.addActionListener(this);
		upperPanel.setLayout(new BoxLayout(upperPanel, BoxLayout.X_AXIS));
		upperPanel.add(new JLabel("Show: "));
		upperPanel.add(selector);
		upperPanel.add(Box.createHorizontalGlue());
		upperPanel
				.add(new JLabel("(C) 2004 ISB. See Help>About for more info"));

		panel.add(upperPanel, BorderLayout.NORTH);
		panel.add(cardPanel, BorderLayout.CENTER);
	}

	public void actionPerformed(ActionEvent ae)
	{
		// selector changed?
		if (ae.getSource() == selector)
		{
			if (selector.getSelectedIndex() == 0) // TIC
			{
				lowerLayout.show(cardPanel, "TIC");
			} else
				lowerLayout.show(cardPanel, "BPI");
		} else if (ae.getSource() == exitMenuItem)
		{
			System.exit(0);
		} else if (ae.getSource() == aboutMenuItem)
		{
			showAboutBox();
		} else if (ae.getSource() == showInfoMenuItem)
		{
			showInfoDialog();
		} else
			System.out.println("Unknown actionsource: " + ae.getSource());
	}

	public void scanClicked(int scanNumber)
	{
		selectScan(scanNumber);
	}

	public void showAboutBox()
	{
		JOptionPane
				.showMessageDialog(
						this,
						"MzXMLViewer "
								+ VERSION_STRING
								+ "\n(C) 2004 ISB (www.systemsbiology.org)\nThis program is FREE software and comes with no warranty.\n"
								+ "This program was made by Mathijs Vogelzang (m_v@dds.nl)\nFor questions about mzXML, contact Patrick Pedrioli (ppatrick@systemsbiology.org)\n",
						"About", JOptionPane.INFORMATION_MESSAGE);

	}

	protected void showInfoDialog()
	{
		new InfoDialog(this, mzXMLFile.getName(), parser.getHeaderInfo());
	}
}