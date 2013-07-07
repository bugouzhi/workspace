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
 * File: * @(#) InfoDialog.java * Author: * Mathijs Vogelzang
 * m_v@dds.nl
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
 * 10-05-2004 Added this header
 * 
 * Created on May 21, 2004
 *  
 ******************************************************************************/
package org.systemsbiology.mzxmlviewer.dialogs;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.table.DefaultTableModel;

import org.systemsbiology.jrap.DataProcessingInfo;
import org.systemsbiology.jrap.MSInstrumentInfo;
import org.systemsbiology.jrap.MZXMLFileInfo;
import org.systemsbiology.jrap.ParentFile;
import org.systemsbiology.jrap.SoftwareInfo;

/**
 * InfoDialog gives information about the header information contained in an 
 * mzXML document.
 * 
 * @author M. Vogelzang
 */
public class InfoDialog extends JDialog
{
	protected JTabbedPane tabbedPane;
	protected MZXMLFileInfo info;

	public InfoDialog(JFrame parent, String fileName, MZXMLFileInfo info)
	{
		super(parent, "Header information for " + fileName, true);

		this.info = info;

		Box mainPanel = new Box(BoxLayout.Y_AXIS);
		mainPanel.setBorder(BorderFactory.createEmptyBorder(11, 11, 11, 11));
		getContentPane().add(mainPanel);

		tabbedPane = new JTabbedPane();
		mainPanel.add(tabbedPane);
		mainPanel.add(Box.createVerticalStrut(11));

		Box buttonPanel = new Box(BoxLayout.X_AXIS);
		buttonPanel.add(Box.createHorizontalGlue());
		JButton ok = new JButton("Ok");
		ok.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent ae)
			{
				dispose();
			}
		});
		buttonPanel.add(ok);
		mainPanel.add(buttonPanel);

		tabbedPane.add(createFileProcPanel(), "Files and Processing");
		tabbedPane.add(createInstrumentPanel(), "MS instrument");

		pack();
		setLocationRelativeTo(parent);
		setVisible(true);
	}

	protected String valueToString(int value)
	{
		switch (value)
		{
			case DataProcessingInfo.YES:
				return "On";
			case DataProcessingInfo.NO:
				return "Off";
			default:
				return "Unknown";
		}
	}

	public JPanel createFileProcPanel()
	{
		JPanel panel = new JPanel(new GridBagLayout());
		panel.setBorder(BorderFactory.createEmptyBorder(11, 11, 11, 11));

		String[] labels = new String[] { "Parent files", "Intensity cutoff",
				"Centroiding", "De-isotoping", "Charge deconvolution",
				"Spot integration", "Sofware used"};

		GridBagConstraints gbc = new GridBagConstraints();
		gbc.insets = new Insets(0, 0, 5, 20);
		gbc.anchor = GridBagConstraints.NORTHWEST;
		for (int i = 0; i < labels.length; i++)
		{
			gbc.gridy = i;
			panel.add(new JLabel(labels[i]), gbc);
		}
		gbc.gridx = 1;
		gbc.gridy = 0;
		gbc.insets = new Insets(0, 0, 0, 0);

		DefaultTableModel parentFileModel = new DefaultTableModel(new String[] {
				"URI", "Sha1-sum", "Type"}, 0);

		JTable parentFileTable = new JTable(parentFileModel);
		parentFileTable.setPreferredScrollableViewportSize(new Dimension(600,
				100));

		ParentFile[] parentFiles = info.getParentFiles();
		for (int i = 0; i < parentFiles.length; i++)
		{
			String type;
			switch (parentFiles[i].getType())
			{
				case ParentFile.TYPE_PROCESSED:
					type = "Processed";
					break;
				case ParentFile.TYPE_RAW:
					type = "Raw";
					break;
				default:
					type = "Unknown";
			}
			parentFileModel.addRow(new String[] { parentFiles[i].getURI(),
					parentFiles[i].getSha1(), type});
		}

		//parentFileTable.doLayout();
		parentFileTable.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);

		panel.add(new JScrollPane(parentFileTable));
		gbc.gridy++;

		DataProcessingInfo dpi = info.getDataProcessing();
		if (dpi == null)
		{
			for (int i = 0; i < labels.length - 1; i++)
			{
				panel.add(new JLabel("Unknown"), gbc);
				gbc.gridy++;
			}
		} else
		{
			DefaultTableModel model = new DefaultTableModel(new String[] {
					"Step", "Software", "Version"}, 0);
			JTable table = new JTable(model);

			SoftwareInfo[] infos = dpi.getSoftwareUsed();
			for (int i = 0; i < infos.length; i++)
			{
				model.addRow(new String[] { infos[i].type, infos[i].name,
						infos[i].version});
			}

			panel.add(new JLabel((dpi.getIntensityCutoff() < 0 ? "Unknown" : ""
					+ dpi.getIntensityCutoff())), gbc);
			gbc.gridy++;
			panel.add(new JLabel(valueToString(dpi.getCentroided())), gbc);
			gbc.gridy++;
			panel.add(new JLabel(valueToString(dpi.getDeisotoped())), gbc);
			gbc.gridy++;
			panel.add(new JLabel(valueToString(dpi.getChargeDeconvoluted())),
					gbc);
			gbc.gridy++;
			panel.add(new JLabel(valueToString(dpi.getSpotIntegration())), gbc);
			gbc.gridy++;
			gbc.fill = GridBagConstraints.BOTH;
			gbc.weightx = 1;
			gbc.weighty = 1;
			table.setPreferredScrollableViewportSize(new Dimension(600, 100));
			panel.add(new JScrollPane(table), gbc);
		}

		return panel;
	}

	public JPanel createInstrumentPanel()
	{
		MSInstrumentInfo mii = info.getInstrumentInfo();

		JPanel panel = new JPanel();
		panel.setBorder(BorderFactory.createEmptyBorder(11, 11, 11, 11));
		if (mii == null)
			panel.add(new JLabel(
					"No information about MS instrument available!"));
		else
		{
			String labels[] = new String[] { "Manufacturer", "Model",
					"Ionisation method", "Mass analyzer", "Detector",
					"Operator name", "Operator URI", "Operator phone",
					"Operator email", "Software name and version"};

			String values[] = new String[] {
					mii.getManufacturer(),
					mii.getModel(),
					mii.getIonization(),
					mii.getMassAnalyzer(),
					mii.getDetector(),
					mii.getOperator() == null ? null
							: mii.getOperator().firstName + " "
									+ mii.getOperator().lastName,
					mii.getOperator() == null ? null : mii.getOperator().URI,
					mii.getOperator() == null ? null
							: mii.getOperator().phoneNumber,
					mii.getOperator() == null ? null : mii.getOperator().email,
					mii.getSoftwareInfo().name + " "
							+ mii.getSoftwareInfo().version};

			panel.setLayout(new GridBagLayout());
			GridBagConstraints gbc = new GridBagConstraints();
			gbc.gridx = gbc.gridy = 0;

			gbc.anchor = GridBagConstraints.NORTHWEST;
			for (int i = 0; i < labels.length; i++)
			{
				gbc.gridx = 0;
				gbc.insets = new Insets(0, 0, 5,20);
				panel.add(new JLabel(labels[i]), gbc);
				gbc.gridx = 1;

				gbc.insets = new Insets(0, 0, 0, 0);
				if (values[i] == null)
					panel.add(new JLabel("Unknown"), gbc);
				else
					panel.add(new JLabel(values[i]), gbc);
				gbc.gridy++;
			}
			gbc.gridx = 2;
			gbc.weightx = 1;
			gbc.weighty = 1;
			panel.add(new JLabel(""), gbc);
		}

		return panel;
	}
}