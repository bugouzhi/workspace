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
 * File: SpectrumComponent.java  
 * Author: Mathijs Vogelzang
 * mvogelzang@systemsbiology.org
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
 * 10-05-2004: adapted this version to JFreechart 0.9.21
 */

package org.systemsbiology.mzxmlviewer.utilities;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.geom.Rectangle2D;
import java.text.NumberFormat;
import java.util.*;

import javax.swing.JPanel;

import org.jfree.chart.*;
import org.jfree.chart.plot.*;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.data.Range;
import org.jfree.data.general.*;
import org.jfree.data.xy.*;
import org.systemsbiology.jrap.*;

/**
 * SpectrumComponent provides a view of a spectrum using JFreeChart.
 * 
 * @author M. Vogelzang
 */
public class SpectrumComponent extends JPanel
{
	static class MyDataSet extends AbstractIntervalXYDataset
	{
		private DatasetGroup group;
		private final static double PEAK_WIDTH = 0.05;
		private Vector listeners;

		float[][] massIntensityList;

		public MyDataSet(float[][] p)
		{
			massIntensityList = p;
			//peaks = new float[][] { { 1, 5}, {4, 3}, {5, 1} };
			listeners = new Vector();
		}

		public void setData(float[][] p)
		{
			massIntensityList = p;
			for (int i = 0; i < listeners.size(); i++)
			{
				((DatasetChangeListener) listeners.get(i))
						.datasetChanged(new DatasetChangeEvent(this, this));
			}
		}

		public Number getStartX(int series, int item)
		{
			return new Double(massIntensityList[0][item] - PEAK_WIDTH / 2);
		}

		public Number getEndX(int series, int item)
		{
			return new Double(massIntensityList[0][item] + PEAK_WIDTH / 2);
		}

		public Number getStartY(int series, int item)
		{
			return new Double(0);
		}

		public Number getEndY(int series, int item)
		{
			return new Double(massIntensityList[1][item]);
		}

		public int getItemCount(int series)
		{
			return massIntensityList[0].length;
		}

		public Number getX(int series, int item)
		{
			return getEndX(series, item);
		}

		public Number getY(int series, int item)
		{
			return getEndY(series, item);
		}

		public int getSeriesCount()
		{
			return 1;
		}

		public String getSeriesName(int series)
		{
			return "Peaks";
		}

		public void addChangeListener(DatasetChangeListener listener)
		{
			listeners.addElement(listener);
		}

		public void removeChangeListener(DatasetChangeListener listener)
		{
			listeners.removeElement(listener);
		}

		public DatasetGroup getGroup()
		{
			return group;
		}

		public void setGroup(DatasetGroup group)
		{
			this.group = group;
		}
	}

	JFreeChart chart;
	MyDataSet data;

	public SpectrumComponent()
	{
		setLayout(new BorderLayout());
	}

	public void setScan(Scan s)
	{
		if (chart == null)
		{
			data = new MyDataSet(s.getMassIntensityList());
			chart = ChartFactory.createXYBarChart("Spectrum @ "
					+ s.getRetentionTime() + " s", "m/z", false, "intensity",
					data, PlotOrientation.VERTICAL, false, false, false);

			//chart.getXYPlot().setDomainAxis(new NumberAxis());
			//chart.getXYPlot().setRangeAxis(new NumberAxis());
			XYBarRenderer renderer = new XYBarRenderer();
			renderer.setSeriesItemLabelsVisible(0, new Boolean(true));
			chart.getXYPlot().setRenderer(renderer);
			removeAll();
			add(new ChartPanel(chart));
			validate();
			repaint();
		}

		NumberFormat format = NumberFormat.getNumberInstance();
		format.setMaximumFractionDigits(1);
		chart.setTitle("Spectrum @ "
				+ format.format(s.getDoubleRetentionTime()) + " s");
		chart.getXYPlot().clearDomainMarkers();
		if (s.getMsLevel() == 2)
			chart.getXYPlot().addDomainMarker(
					new ValueMarker(s.getPrecursorMz()));
		data.setData(s.getMassIntensityList());
	}

	public Dimension getPreferredSize()
	{
		return new Dimension(640, 480);
	}

	public Dimension getMinimumSize()
	{
		return getPreferredSize();
	}

	public interface ClickListener
	{
		public void scanClicked(int scanNumber);
	}

	private static class CMListAdapter implements ChartMouseListener
	{
		private ClickListener cl;
		private double[] retTimes;

		CMListAdapter(ClickListener cl, ScanHeader[] headers)
		{
			this.cl = cl;
			retTimes = new double[headers.length];
			for (int i = 0; i < retTimes.length; i++)
				retTimes[i] = headers[i].getDoubleRetentionTime();
		}

		public void chartMouseClicked(ChartMouseEvent event)
		{
			JFreeChart chart = event.getChart();
			XYPlot plot = chart.getXYPlot();
			MouseEvent mouseEvent = event.getTrigger();
			ChartPanel chartPanel = (ChartPanel) mouseEvent.getSource();

			Rectangle2D rect = chartPanel.getScaledDataArea();

			int x = mouseEvent.getX();
			if (rect.getMinX() <= x && rect.getMaxX() >= x)
			{
				Range range = plot.getDomainAxis().getRange();
				double rangeRatio = (x - rect.getMinX()) / (rect.getWidth());
				int number;
				double time = range.getLowerBound() + rangeRatio
						* range.getLength();
				if (time <= retTimes[0])
					number = 1;
				else if (time >= retTimes[retTimes.length - 1])
					number = retTimes.length;
				else
				{
					number = Arrays.binarySearch(retTimes, time);
					if (number < 0)
						number = -(number + 1) + 1;
					else
						number++;
				}
				cl.scanClicked(number);
			}
		}

		public void chartMouseMoved(ChartMouseEvent event)
		{
			// Do nothing
		}
	}

	public static Component getBPIComponent(ScanHeader[] scanHeaders,
			ClickListener cl)
	{
		XYSeries dataSeries = new XYSeries("TIC");
		for (int i = 0; i < scanHeaders.length; i++)
			dataSeries.add(scanHeaders[i].getDoubleRetentionTime(),
					scanHeaders[i].getTotIonCurrent());
		XYSeriesCollection dataSet = new XYSeriesCollection(dataSeries);
		JFreeChart BPIChart = ChartFactory.createXYAreaChart(null,
				"Retention time (s)", "Intensity", dataSet,
				PlotOrientation.VERTICAL, false, false, false);
		ChartPanel cPanel = new ChartPanel(BPIChart);
		cPanel.setPreferredSize(new Dimension(640, 300));
		cPanel.addChartMouseListener(new CMListAdapter(cl, scanHeaders));
		return cPanel;
	}

	public static Component getTICComponent(ScanHeader[] scanHeaders,
			ClickListener cl)
	{
		XYSeries dataSeries = new XYSeries("BPI");
		for (int i = 0; i < scanHeaders.length; i++)
			dataSeries.add(scanHeaders[i].getDoubleRetentionTime(),
					scanHeaders[i].getBasePeakIntensity());
		XYSeriesCollection dataSet = new XYSeriesCollection(dataSeries);
		JFreeChart BPIChart = ChartFactory.createXYAreaChart(null,
				"Retention time (s)", "Intensity", dataSet,
				PlotOrientation.VERTICAL, false, false, false);

		ChartPanel cPanel = new ChartPanel(BPIChart);
		cPanel.setPreferredSize(new Dimension(640, 300));
		cPanel.addChartMouseListener(new CMListAdapter(cl, scanHeaders));
		return cPanel;

	}
}