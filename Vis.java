import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import javax.swing.*;
import javax.imageio.*;
import java.awt.image.*;
import javax.swing.filechooser.*;
import java.util.*;
import java.io.*;

public class Vis extends JApplet implements ActionListener{
	public static String version="8th Feb 2016";
	final static Color bg = Color.white;
	final static Color fg = Color.black;
	Dimension totalSize;
	public BufferedImage img;
	public int size=5000;
	int numParticles=0;
	int[][] OriginalMatrix;
	int mag=1; // Magnification factor
	public CA ca;
	int counter=0; // Timesteps counter
	int defaultView=0; // Cells view
	int timestep=0;
	boolean movie=true;
	public boolean play=true;
	JLabel label;
	String rootDir=".";


	public Vis (JLabel l, float densityV)
	{
		label = l;
		ca = new CA();
		size=ca.size;
		mag = ca.mag;
		img = new BufferedImage (ca.size*ca.mag,1*ca.size*ca.mag,BufferedImage.TYPE_INT_ARGB);
	}

	public final BufferedImage scale(double scale, BufferedImage srcImg)
	{
		if (scale == 1)  return srcImg;
   		AffineTransformOp op = new AffineTransformOp( AffineTransform.getScaleInstance( scale, scale), null);
      return op.filter(srcImg, null);
	}

	public void paint(Graphics g) {
      Graphics2D g2 = (Graphics2D) g;
			BufferedImage img2=scale (mag,img);
			g2.drawImage(img2,null,null);
			if (ca.size!=size) {
				img = new BufferedImage (ca.size*ca.mag,ca.size*ca.mag,BufferedImage.TYPE_INT_ARGB);
				ca.size=size;
			}
	}

	public void nextTimeStep ()
	{
		ca.nextTimeStep();
		paintAll();
		timestep++;

        if (timestep==1000) ca.chemo=true;
	}


	public void paintAll()
	{
		int [][] lattice = ca.getCells();

		BufferedImage all = new BufferedImage (ca.size*2,ca.size*1,BufferedImage.TYPE_INT_ARGB);
		BufferedImage c= getCells();
		for (int i=0;i<ca.size;i++)
			for (int j=0;j<size;j++) {
				all.setRGB(i,j,c.getRGB(i,j));
				img.setRGB(i,j,c.getRGB(i,j));
			}

		repaint();
		if (movie) try {
			File dir = new File (rootDir+"/images");
			dir.mkdir ();
			File fileAll = new File (rootDir+"/images/all"+counter+".png");
			ImageIO.write(img,"png",fileAll);
		} catch (Exception e) {
			e.printStackTrace();
		}

		counter++;
	}

	public BufferedImage getCells()
	{
		int [][] lattice = ca.getCells();
		int [][] vasculature = ca.getVasculature();
		BufferedImage result = new BufferedImage (ca.size,ca.size,BufferedImage.TYPE_INT_ARGB);

		for (int i=0;i<ca.size;i++)
			for (int j=0;j<ca.size;j++) {
				int val=0;
				if (lattice[i][j]==0) val=Color.white.getRGB();
				else if (lattice[i][j]<=10) val=Color.black.getRGB();
				else if (lattice[i][j]<=20) val=Color.darkGray.getRGB();
				else if (lattice[i][j]<=30) val=Color.gray.getRGB();
                else if (lattice[i][j]<=40) val=Color.lightGray.getRGB();
                else if (lattice[i][j]<=50) val=Color.cyan.getRGB();
				else if (lattice[i][j]<=60) val=Color.yellow.getRGB();
				else if (lattice[i][j]<=70) val=Color.orange.getRGB();
				else if (lattice[i][j]<=80) val=Color.magenta.getRGB();
				else if (lattice[i][j]>80) val=Color.red.getRGB();
				result.setRGB(i,j,val);
			}
		return result;
	}

	public void actionPerformed (ActionEvent e) {
		String res = e.getActionCommand();
		repaint();
	}
	public void mouseExited (MouseEvent e) {}
	public void mouseEntered (MouseEvent e) {}
	public void mouseReleased (MouseEvent e) {}
	public void mousePressed (MouseEvent e) {}


	public static void main(String args[]) {
		int mag=1;
		int maxTS=1000;
		boolean movie;

		System.err.println ("# Vis version:"+Vis.version);
		System.err.println ("# CA version:"+CA.version);

		float densityV=0.04f;
		if (args.length==2) {
			maxTS=Integer.parseInt (args[0]);
			densityV=Float.parseFloat(args[1]);
		} else {
			System.err.println ("Arguments needed: timesteps, densityV");
			System.exit(-1);
		}
		// Let's deal now with the main window
		JFrame f = new JFrame("Chemotherapy CA Visualisation "+CA.version);
	        f.addWindowListener(new WindowAdapter() {
	        	public void windowClosing(WindowEvent e) {System.exit(0);}
	        });
		JLabel label = new JLabel ("Timesteps");


		Vis m = new Vis (label,densityV);
		JPanel buttonPanel = new JPanel(); //use FlowLayout

		f.getContentPane().add("Center", m);
		f.getContentPane().add("North",label);
		f.pack();
		f.setSize(new Dimension(m.ca.mag*m.size,m.ca.mag*m.size+80));
		f.show();

		boolean finished=false;
		int ts=0;
		while (ts!=maxTS) {
			if (m.play) {
				m.nextTimeStep();
				label.setText ("Timestep: "+ts);
				ts++;
			}
		}
		System.exit(-1);
	}

}
