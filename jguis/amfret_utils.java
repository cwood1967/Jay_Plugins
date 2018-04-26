package jguis;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.frame.RoiManager;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;

import java.io.File;

import jalgs.algutils;
import jalgs.gui_interface;
import jalgs.jdataio;
import jalgs.jstatistics;
import jalgs.jfit.fit_gaussian;
import jalgs.jfit.fit_stretched_exp;

public class amfret_utils implements gui_interface{

	public static void main(String[] args){
		//here we batch process a folder outside of imagej
		//the question is how to pass all of the parameters--all as args?
		//alternatively, can fix all to defaults and change the code if necessary
		//compromise: pass the conc parameters and fix everything else
		if(args[0].equals("--getgate")) {
			(new amfret_utils()).getGate(args[1],args[2],args[3],args[4],Float.parseFloat(args[5]),Float.parseFloat(args[6]),-0.2f,1.0f,1.5f,99.0f,99.0f,6,false);
			return;
		}
		//System.out.println(table_tools.print_string_array(args,1));
		String[] fcsfiles=(new jdataio()).get_sorted_string_list(args[0],null);
		String[] col_labels={"title","acceptor","c^2","Iter","baseline","amp","EC50","alpha","xshift","EC50_errs","alpha_errs","totcells","fretcells","bimodal_metric","f_gate","delta","delta_errs"};
		System.out.println(table_tools.print_string_array(col_labels,1));
		amfret_utils custom=new amfret_utils();
		float minamfret=-0.2f; float maxamfret=1.0f; int mincells=1; float mincrop=1.5f;
		if(args.length>5) {
			//passing other variables
			minamfret=Float.parseFloat(args[5]);
			if(args.length>6) maxamfret=Float.parseFloat(args[6]);
			if(args.length>7) mincells=(int)Float.parseFloat(args[7]);
			if(args.length>8) mincrop=Float.parseFloat(args[8]);
		}
		try{
			BufferedWriter b=new BufferedWriter(new FileWriter(new File(args[1]+"Stretched Exp Fits.xls")));
			b.write(table_tools.print_string_array(col_labels,0)+"\n");
			for(int i=0;i<fcsfiles.length;i++){
				if(!fcsfiles[i].endsWith(".fcs")) continue;
				//assuming that indir and outdir are the same
				List<String> output=custom.exec(args[0],fcsfiles[i],args[1],args[2],Float.parseFloat(args[3]),Float.parseFloat(args[4]),minamfret,maxamfret,mincells,mincrop,0.02f,0.95f,false);
				if(output!=null){
					System.out.println(table_tools.print_string_array(output,1));
					b.write(table_tools.print_string_array(output,0));
				}
			}
			b.close();
		} catch (IOException e){
			System.out.println("error writing file");
		}
		return;
	}

	public List<String> exec(String directory,String name,String outdir,String roipath,float minconc,float maxconc,float minamfret,float maxamfret,
		int mincells,float startcrop,float minbimodefrac,float maxbimodefrac,boolean showplots){
		//get the gate roi
		boolean loadroi=true;
		if(roipath==null) loadroi=false;
		Roi[] rois=null;
		if(loadroi){
			rois=new Roi[1];
			rois[0]=multi_roi_writer.readRoi(roipath);
		} else {
			RoiManager rman=RoiManager.getInstance();
			if(rman==null){
				IJ.error("Need gate roi in roi manager");
				return null;
			}
			rois=rman.getRoisAsArray();
		}
		//now load the fcs data
		Object[] data=(new import_flowcyte()).getFCSFile(directory+name);
		//now calculate the acceptor and amfret data sets
		//find the desired column names
		String shortname=name.substring(0,name.length()-4);
		float[][] rowdata=(float[][])data[1];
		if(rowdata==null){
			if(showplots) IJ.log("No data in the fcs file");
			List<String> blankout=new ArrayList<String>();
			blankout.add(shortname);
			for(int i=1;i<17;i++) blankout.add("NaN");
			return blankout;
		}
		if(rowdata.length<mincells){
			if(showplots) IJ.log("Not enough cells in the fcs file");
			List<String> blankout=new ArrayList<String>();
			blankout.add(shortname);
			for(int i=1;i<17;i++) blankout.add("NaN");
			return blankout;
		}
		int ncells=rowdata.length;
		String[] colnames=(String[])data[0];
		String accname="FL10-A"; int acccol=0;
		String fretname="FL04-A"; int fretcol=0;
		for(int i=0;i<colnames.length;i++){
			if(accname.equals(colnames[i])){
				acccol=i;
				break;
			}
		}
		for(int i=0;i<colnames.length;i++){
			if(fretname.equals(colnames[i])){
				fretcol=i;
				break;
			}
		}
		float[] acceptor=new float[rowdata.length];
		float[] amfret=new float[rowdata.length];
		for(int i=0;i<rowdata.length;i++){
			acceptor[i]=rowdata[i][acccol]/1000000.0f;
			amfret[i]=rowdata[i][fretcol]/rowdata[i][acccol];
		}
		float acceptoravg=jstatistics.getstatistic("Avg",acceptor,null);

		//now create the 2D AmFRET histogram
		Plot2DHist amfrethist=new Plot2DHist("Acceptor","AmFRET",acceptor,amfret,null);
		amfrethist.setLogAxes(true,false);
		float[] templims=amfrethist.getLimits();
		amfrethist.setLimits(minconc,maxconc,minamfret,maxamfret,templims[4],templims[5]);
		amfrethist.setBinSize(4);
		amfrethist.intautoscale();
		String histname=name.substring(0,name.length()-4)+"_hist.pw2";
		if(showplots) (new PlotWindow2DHist(histname,amfrethist)).draw();

		//now get the acceptor histogram 99.9th percentile
		float uplim=jstatistics.getstatistic("Percentile",acceptor,new float[]{99.9f});
		if(showplots) IJ.log("fitting upper limit = "+uplim);
		
		//now get the grayscale AmFRET image and gate it
		FloatProcessor histimage=amfrethist.getHistImage();
		float[] histpix=(float[])histimage.getPixels();
		int histwidth=histimage.getWidth();
		int histheight=histimage.getHeight();
		float[] totpop=new float[histwidth];
		float[] fretpop=new float[histwidth];
		float[] fracfret=new float[histwidth];
		float totcells=0.0f;
		float fretcells=0.0f;
		for(int i=0;i<histwidth;i++){
			for(int j=0;j<histheight;j++){
				totpop[i]+=histpix[i+j*histwidth];
				if(!rois[0].contains(i,j)) fretpop[i]+=histpix[i+j*histwidth];
			}
			totcells+=totpop[i]/16.0f; //divide by 16 because bins are 4 x 4
			fretcells+=fretpop[i]/16.0f;
			if(totpop[i]>0.0f) fracfret[i]=fretpop[i]/totpop[i];
		}
		float fretfrac=fretcells/totcells;

		//now calculate the x axis
		float logxmin=(float)Math.log(minconc);
		float logxmax=(float)Math.log(maxconc);
		float logxinc=(logxmax-logxmin)/(float)(histwidth-1);
		float[] xvals=new float[histwidth];
		for(int i=0;i<histwidth;i++){
			float logval=logxmin+logxinc*(float)i;
			xvals[i]=(float)Math.exp(logval);
		}
		Plot4 gatefracplot=new Plot4("Acceptor","Fraction Assembled",xvals,fracfret);
		templims=gatefracplot.getLimits();
		gatefracplot.setLimits(templims[0],templims[1],0.0,1.0);
		gatefracplot.setLogAxes(true,false);
		String gatefracplotname=name.substring(0,name.length()-4)+"_gatefraction.pw2";
		if(showplots) (new PlotWindow4(gatefracplotname,gatefracplot)).draw();

		//now crop the fracfret profile
		int counter1=-1;
		do{
			counter1++;
		}while(xvals[counter1]<startcrop);
		int counter2=counter1;
		do{
			counter2++;
		}while(counter2<(xvals.length-1) && xvals[counter2]<uplim);
		int newlength=counter2-counter1+1;
		float[] newxvals=(float[])algutils.get_subarray(xvals,counter1,newlength);
		float[] newfracfret=(float[])algutils.get_subarray(fracfret,counter1,newlength);
		
		//now fit to a stretched exponential
		fit_stretched_exp fitclass=new fit_stretched_exp(newxvals,newfracfret,startcrop,uplim,this);
		if(showplots) IJ.log("Starting Parameters");
		if(showplots) IJ.log(""+fitclass.sparams[0]+" , "+fitclass.sparams[1]+" , "+fitclass.sparams[2]+" , "+fitclass.sparams[3]+" , "+fitclass.sparams[4]);
		int[] fixes={1,1,0,0,1};
		Object[] fitresults=fitclass.autoFit(fixes,fitclass.sparams,null,true,false);
		float[] fit=(float[])fitresults[2];
		double[] params=(double[])fitresults[0]; //params are baseline, amp, ec50, alpha, and xcmin
		double[] stats=(double[])fitresults[1];
		float[] errs=(float[])fitresults[3];
		if(showplots) IJ.log("Fit Results");
		if(showplots) IJ.log(""+params[0]+" , "+params[1]+" , "+params[2]+" , "+params[3]+" , "+params[4]);
		if(showplots) IJ.log("Fit Errs (unfixed params and c2)");
		if(showplots) IJ.log(""+errs[0]+" , "+errs[1]+" , "+errs[2]); //this is ec50, alpha, and chi squared
		Plot4 fitplot=new Plot4("Acceptor","Fraction Assembled",new float[][]{newxvals,newxvals},new float[][]{newfracfret,fit},null);
		templims=fitplot.getLimits();
		fitplot.setLimits(templims[0],templims[1],0.0,1.0);
		fitplot.setLogAxes(true,false);
		String fitname=name.substring(0,name.length()-4)+"_gatefraction_crop.pw2";
		if(showplots) (new PlotWindow4(fitname,fitplot)).draw();

		//if the fretting fraction is between the preset limits, check for bimodality
		float bimodemetric=0.0f;
		if(fretfrac>=minbimodefrac && fretfrac<=maxbimodefrac){
			float transwidth=(float)(params[2]/(params[3]*0.3466));
			float transwinstart=(float)(params[2]-0.25*transwidth);
			float transwinend=(float)(params[2]+0.25*transwidth);
			if(transwinstart<startcrop) transwinstart=startcrop;
			bimodemetric=detect2DHistBiphasic(acceptor,amfret,transwinstart,transwinend,minamfret,maxamfret,1.0f,false);
		}

		//finally append all parameters to a table and save the histograms and plots
		double delta=1.0/params[3];
		double deltaerr=delta*delta*errs[2];

		List<String> output=new ArrayList<String>();
		output.add(shortname);
		output.add(""+acceptoravg);
		output.add(""+(float)stats[1]);
		output.add(""+(int)stats[0]);
		for(int i=0;i<params.length;i++) output.add(""+(float)params[i]);
		output.add(""+(float)errs[0]);
		output.add(""+(float)errs[1]);
		output.add(""+totcells);
		output.add(""+fretcells);
		output.add(""+bimodemetric);
		output.add(""+fretfrac);
		output.add(""+(float)delta);
		output.add(""+(float)deltaerr);
		amfrethist.saveplot2file(outdir+histname);
		gatefracplot.saveplot2file(outdir+gatefracplotname);
		fitplot.saveplot2file(outdir+fitname);
		//make a side-by-side merge of the amfret and gate images as an overview
		ColorProcessor amfretcp=amfrethist.getProcessor();
		ColorProcessor fitcp=fitplot.getProcessor();
		ColorProcessor overviewcp=mergeImages(amfretcp,fitcp);
		ImagePlus overviewimp=new ImagePlus("Overview",overviewcp);
		String overviewname=name.substring(0,name.length()-4)+"_overview.png";
		IJ.save(overviewimp,outdir+overviewname);
		return output;
	}

	public ColorProcessor mergeImages(ColorProcessor cp1,ColorProcessor cp2){
		int width1=cp1.getWidth(); int height1=cp1.getHeight(); int[] pix1=(int[])cp1.getPixels();
		int width2=cp2.getWidth(); int height2=cp2.getHeight(); int[] pix2=(int[])cp2.getPixels();
		int newwidth=width1+width2;
		int newheight=(int)Math.max(height1,height2);
		int[] newpix=new int[newwidth*newheight];
		for(int i=0;i<newwidth*newheight;i++) newpix[i]=0xffffffff;
		for(int i=0;i<height1;i++){
			System.arraycopy(pix1,i*width1,newpix,i*newwidth,width1);
		}
		for(int i=0;i<height2;i++){
			System.arraycopy(pix2,i*width2,newpix,i*newwidth+width1,width2);
		}
		return (new ColorProcessor(newwidth,newheight,newpix));
	}

	public float detect2DHistBiphasic(float[] acceptor,float[] amfret,float startx,float endx,float starty,float endy,float peakwindow,boolean showfit){
		float[] newyvals=new float[amfret.length];
		int nsel=0;
		for(int i=0;i<amfret.length;i++){
			if(acceptor[i]>=startx && acceptor[i]<=endx && amfret[i]>=starty && amfret[i]<=endy){
				newyvals[nsel]=amfret[i];
				nsel++;
			}
		}
		float[] newyvals2=(float[])algutils.get_subarray(newyvals,0,nsel);
		int nbins=64;
		float bininc=(endy-starty)/(float)nbins;
		float[] hist=new float[nbins];
		float[] histxvals=new float[nbins];
		for(int i=0;i<newyvals2.length;i++){
			int temp=(int)Math.floor((newyvals2[i]-starty)/bininc);
			hist[temp]++;
		}
		for(int i=0;i<nbins;i++){
			histxvals[i]=starty+bininc*(float)i;
		}
		//now fit the histogram to a 1D gaussian
		fit_gaussian fg=new fit_gaussian();
		double[] params=new double[4];
		int[] fixes={0,0,0,0};
		double c2=0.0f;
		int iterations=0;
		double[] stats=new double[2];
		float[] fit=fg.run1DFit(histxvals,hist,params,stats,null,fixes);
		if(showfit){
			new PlotWindow4("Y Hist Fit","y","freq",new float[][]{histxvals,histxvals},new float[][]{hist,fit},null).draw();
		}
		float[] resid=new float[hist.length];
		for(int i=0;i<hist.length;i++){
			resid[i]=hist[i]-fit[i];
		}
		//params are baseline, xc, stdev, amp
		float maxresid=0.0f;
		for(int i=0;i<hist.length;i++){
			if(Math.abs(histxvals[i]-(float)params[1])>(peakwindow*(float)params[2])){
				if(resid[i]>maxresid) maxresid=resid[i];
			}
		}
		maxresid/=(float)params[3];
		return maxresid;
	}
	
	public void getGate(String directory,String name,String outdir,String roiname,float minconc,float maxconc,float minamfret,float maxamfret,
			float startcrop,float uppercentile,float gateper,int gateshift,boolean showplots) {
		//now load the fcs data
		Object[] data=(new import_flowcyte()).getFCSFile(directory+name);
		//now calculate the acceptor and amfret data sets
		//find the desired column names
		float[][] rowdata=(float[][])data[1];
		if(rowdata==null){
			IJ.log("No data in the fcs file");
			return;
		}
		int ncells=rowdata.length;
		String[] colnames=(String[])data[0];
		String accname="FL10-A"; int acccol=0;
		String fretname="FL04-A"; int fretcol=0;
		for(int i=0;i<colnames.length;i++){
			if(accname.equals(colnames[i])){
				acccol=i;
				break;
			}
		}
		for(int i=0;i<colnames.length;i++){
			if(fretname.equals(colnames[i])){
				fretcol=i;
				break;
			}
		}
		float[] acceptor=new float[rowdata.length];
		float[] amfret=new float[rowdata.length];
		for(int i=0;i<rowdata.length;i++){
			acceptor[i]=rowdata[i][acccol]/1000000.0f;
			amfret[i]=rowdata[i][fretcol]/rowdata[i][acccol];
		}

		//now create the 2D AmFRET histogram
		Plot2DHist amfrethist=new Plot2DHist("Acceptor","AmFRET",acceptor,amfret,null);
		amfrethist.setLogAxes(true,false);
		float[] templims=amfrethist.getLimits();
		amfrethist.setLimits(minconc,maxconc,minamfret,maxamfret,templims[4],templims[5]);
		amfrethist.setBinSize(4);
		amfrethist.intautoscale();
		String histname=name.substring(0,name.length()-4)+"_hist.pw2";
		if(showplots) (new PlotWindow2DHist(histname,amfrethist)).draw();

		//now get the acceptor histogram upper gate percentile
		float uplim=jstatistics.getstatistic("Percentile",acceptor,new float[]{uppercentile});
		IJ.log("gating upper limit = "+uplim);
		
		//now get the grayscale AmFRET image and determine the gate from it
		FloatProcessor histimage=amfrethist.getHistImage();
		int width=histimage.getWidth(); int height=histimage.getHeight();
		float[] pix=(float[])histimage.getPixels();
		float xmin=minconc; float xmax=maxconc; float ymin=minamfret; float ymax=maxamfret;
		float lowlim=startcrop;

		//generate the x axis (logarithmic)
		float logxmin=(float)Math.log(xmin);
		float logxmax=(float)Math.log(xmax);
		float logxinc=(logxmax-logxmin)/(float)(width-1);
		float[] xvals=new float[width];
		for(int i=0;i<width;i++){
			float logval=logxmin+logxinc*(float)i;
			xvals[i]=(float)Math.exp(logval);
		}
		//now bin by 4
		float[] xvalsbinned=new float[width/4];
		for(int i=0;i<width;i+=4){
			for(int j=0;j<4;j++){
				xvalsbinned[i/4]+=xvals[i+j];
			}
			xvalsbinned[i/4]*=0.25f;
		}
		//now get the y vals
		float[] yvals=new float[height];
		float yinc=(ymax-ymin)/(float)(height-1);
		for(int i=0;i<height;i++){
			//yvals[i]=ymin+yinc*(float)i;
			yvals[i]=(float)i;
		}
		//and bin them by 4
		float[] yvalsbinned=new float[height/4];
		for(int i=0;i<height;i+=4){
			for(int j=0;j<4;j++){
				yvalsbinned[i/4]+=yvals[i+j];
			}
			yvalsbinned[i/4]*=0.25f;
		}
		//now start on the first full bin (first above lowlim) and end on the last full bin (last-1 above uplim)
		int xpos=0;
		while(xvalsbinned[xpos]<lowlim){
			xpos++;
		}
		int firstbin=xpos;
		while(xpos<xvalsbinned.length && xvalsbinned[xpos]<uplim){
			xpos++;
		}
		int lastbin=xpos-1;
		//now go through and get the gate value for each bin
		float[][] hists=new float[width/4][height/4];
		float[][] dupyvals=new float[width/4][height/4];
		for(int i=0;i<width;i+=4){
			dupyvals[i/4]=yvalsbinned;
			for(int j=(height-1);j>=0;j-=4){
				hists[i/4][(height-j+1)/4]=pix[i+j*width];
			}
		}
		//optionally plot the histograms
		//new PlotWindow4("Histograms","bin","hist",dupyvals,hists,null).draw();
		//and finally find the appropriate percentile for each histogram between firstbin and lastbin
		float[] gatevals=new float[width/4];
		for(int i=firstbin;i<=lastbin;i++){
			gatevals[i]=getHistPercentile(hists[i],yvalsbinned,gateper);
		}
		//now fill in the lower and upper bins
		for(int i=0;i<firstbin;i++) gatevals[i]=gatevals[firstbin];
		for(int i=(lastbin+1);i<gatevals.length;i++) gatevals[i]=gatevals[lastbin];
		//should we smooth this?--yes
		smooth(gatevals,5);
		//now plot the gate profile for reference
		if(showplots) new PlotWindow4("Gate Profile","x","y",xvalsbinned,gatevals).draw();
		//now create the roi
		//need to map the gate values back to original bins
		int nroipoints=width/4+4;
		int[][] coords=new int[2][nroipoints];
		for(int i=0;i<width/4;i++){
			coords[0][i+4]=i*4+2;
			//coords[1][i+4]=height-2-(int)((gatevals[i]-ymin)/yinc);
			coords[1][i+4]=height-7-(int)gatevals[i]-gateshift;
		}
		coords[0][0]=width-1; coords[1][0]=coords[1][nroipoints-1];
		coords[0][1]=width-1; coords[1][1]=height-1;
		coords[0][2]=0; coords[1][2]=height-1;
		coords[0][3]=0; coords[1][3]=coords[1][4];
		Roi roi=new PolygonRoi(coords[0],coords[1],nroipoints,Roi.POLYGON);
		if(showplots){
			RoiManager rman=RoiManager.getInstance();
			if(rman==null) rman=new RoiManager();
			rman.addRoi(roi);
		}
		multi_roi_writer.writeRoi(roi,outdir+roiname);
	}
	
	public void smooth(float[] traj,int width){
		//does a y smooth on the trajectory
		//width should be an odd number
		float[] traj2=traj.clone();
		int shift=width/2;
		for(int i=shift;i<(traj.length-shift-1);i++){
			traj2[i]=0.0f;
			for(int j=(i-shift);j<(i-shift+width);j++){
				traj2[i]+=traj[j];
			}
			traj2[i]/=(float)width;
		}
		for(int i=0;i<traj.length;i++) traj[i]=traj2[i];
	}

	public float getHistPercentile(float[] hist,float[] axis,float percentile){
		//this returns a hist percentile axis units
		//note that percentile is a percentile, not a fraction
		//get the normalized cumulative histogram
		float sum=jstatistics.getstatistic("Sum",hist,null);
		float[] cum=new float[hist.length];
		cum[0]=hist[0]/sum;
		for(int i=1;i<hist.length;i++){
			cum[i]=cum[i-1]+hist[i]/sum;
		}
		float frac=percentile/100.0f;
		//finally get the percentile
		for(int i=0;i<hist.length;i++){
			if(cum[i]==frac){
				return axis[i];
			} else if(cum[i]>frac){
				float rem=cum[i]-frac;
				float interp=axis[i-1]+rem*(axis[i]-axis[i-1]);
				return interp;
			}
		}
		return axis[hist.length-1];
	}

	public void showMessage(String message){
		IJ.log(message);
	}

	public void showProgress(int currpos,int finalpos){
		IJ.showProgress(currpos,finalpos);
	}

}
