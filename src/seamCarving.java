
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.awt.image.BufferedImage;

import javax.imageio.ImageIO;

public class seamCarving {
	public static int Pcounter = 0;//number of calculation spared by memoization of p-values
	public static int Gcounter = 0;//number of calculation spared by memoization of greyscale
	
	public static void main(String args[]) {
		forwardEnergyTest();
		
	}
	public static void forwardEnergyTest(){
		BufferedImage img = readImage("images\\bikerider.jpg");
		for(int i=0;i<300;i++){
			System.out.println(i);
			double[][][] costs = forwardEnergyCost(img);
			double[][] dynamicForwardMap = dynamicForwardEnergyMap(costs);
			double[][] Emap = createEnegryMap(img, img.getWidth(),img.getHeight());
			double[][] dynamicEnergyMap = dynamicMap(Emap,0);
			double [][]Dmap = weightedAverageOfMatrices(dynamicForwardMap,dynamicEnergyMap,0.5);
			int[] seam =  chooseSeam(Dmap,0);
			img= deleteSeam( img, 0, seam);
		
		}
			
		
		File f = new File("images\\result.jpg");
		try{
	    ImageIO.write(img, "jpg", f);}
		catch(IOException e){
		      System.out.println("Error: "+e);}
	}
	public static void entropyTest(){
		BufferedImage image = null;
		image = readImage("images\\hen.jpg");
		System.out.println(pixelEntropy(image,20,7,null,null,0));
		/*double[][] Gmap = imageToGreyScaleMatrix(image);
		double[][] Entmap = returnEntropy(image,0);
		System.out.println(Gmap.length+","+Gmap[0].length);
		exportMatrix(Gmap,"Gmap.txt");
		exportMatrix(Entmap,"Entmap.txt");*/
		
	}
	public static void entropyTimeTest(){
		 File f=new File("images\\E.jpg");
			BufferedImage img=null;
		double[] times = new double[10];
			try {
				 img = ImageIO.read(f);
			} catch (IOException e) {
				
				e.printStackTrace();
			}
			
			int width = img.getWidth();
			int height = img.getHeight();
		
				times[0] = System.nanoTime()/1000000000;
				double energy_map[][] = createEnegryMap(img, width,height);
				//System.out.println("finish creating energy map");
				times[1] = System.nanoTime()/1000000000;
				System.out.println("energy time:"+(times[1]-times[0]));
				ExportMatrixAsImage(energy_map,img,"images\\energyMap.jpg");
				
				times[2] = System.nanoTime()/1000000000;
				double[][] energy_map_after_entropy =addEntropy(energy_map,img,1);
				times[3] = System.nanoTime()/1000000000;
				System.out.println("entorpy time:"+(times[3]-times[2]));
				ExportMatrixAsImage(energy_map_after_entropy,img,"images\\energyMapWithEntropy.jpg");
				System.out.println("memoization Ind:"+1);
				System.out.println("Pcounter:"+Pcounter);
				System.out.println("Gcounter:"+Gcounter+"\n");
				
				Pcounter=Gcounter=0;
				times[2] = System.nanoTime()/1000000000;
				energy_map_after_entropy =addEntropy(energy_map,img,3);
				times[3] = System.nanoTime()/1000000000;
				System.out.println("entorpy time:"+(times[3]-times[2]));
				ExportMatrixAsImage(energy_map_after_entropy,img,"images\\energyMapWithEntropy.jpg");
				System.out.println("memoization Ind:"+3);
				System.out.println("Pcounter:"+Pcounter);
				System.out.println("Gcounter:"+Gcounter+"\n");
	}

	
	public static void chooseAndDeleteSeamTest(){
		BufferedImage img = readImage("images\\catcat.jpg");
		for(int i=0;i<500;i++){
			System.out.println(i);
			double[][] Emap = createEnegryMap(img, img.getWidth(),img.getHeight());
			double[][] Dmap = dynamicMap(Emap,0);
			int[] seam =  chooseSeam(Dmap,0);
			img= deleteSeam( img, 0, seam);
			}
		
		File f = new File("images\\result.jpg");
		try{
	    ImageIO.write(img, "jpg", f);}
		catch(IOException e){
		      System.out.println("Error: "+e);
		    }
		
		
	}
	public static double[][] weightedAverageOfMatrices(double[][] m1,double[][]m2,double w){
		/**0<=w<=1 s.t if w=1 then res=m1**/
		if(m1.length!=m2.length||m1[0].length!=m2[0].length){
			System.out.println("Error: unmatched sizes!");
			return null;
		}
		int n= m1.length;
		int m = m1[0].length;
		double[][] res = new double[n][m];
		for(int i=0;i<n;i++){
			for(int j=0;j<m;j++){
				res[i][j]=w*m1[i][j]+(1-w)*m2[i][j];
			}
		}
		return res;
		
	}
	public static double[][][] forwardEnergyCost(BufferedImage img){
		/**res[i][j][0]=CL, res[i][j][1]=CU,res[i][j][2]=CR**/
		int w = img.getHeight();
		int h = img.getWidth();
		double [][][] res = new double[h][w][3];
		int i= 0;
		int j=0;
		//Upper edges has no forward cost
		res[0][0][0]=0;
		res[0][w-1][0]=0;
		//handling the first row
		for(j=1;j<w-1;j++){
			res[i][j][0]=Math.abs(greyscale(img,i,j-1)-greyscale(img,i,j-1));

			
		}

		//handling the first column (no Cost-left)
		j=0;
		for(i=1;i<h;i++){

			res[i][j][1]=0;
			res[i][j][2]=Math.abs(greyscale(img,i-1,j)-greyscale(img,i,j+1));
		}
		//handling the last column (no Cost-right)
		j=w-1;
		for(i=1;i<h;i++){
			res[i][j][1]=0;
			res[i][j][0]=Math.abs(greyscale(img,i-1,j)-greyscale(img,i,j-1));
		}
		// handling the rest of the matrix
		for(i=1;i<h;i++){
			for(j=1;j<w-1;j++){
				res[i][j][0]=c(img,i,j,'L');
				res[i][j][1]=c(img,i,j,'U');
				res[i][j][2]=c(img,i,j,'R');
			}
		}
		return res;
	}
	public static double c(BufferedImage img,int i,int j,char direction){
		if(direction=='L'){
			return Math.abs(greyscale(img,i,j+1)-greyscale(img,i,j-1))+Math.abs(greyscale(img,i-1,j)-greyscale(img,i,j-1));
		}
		if(direction=='U'){
			return Math.abs(greyscale(img,i,j+1)-greyscale(img,i,j-1));
		}
		if(direction=='R'){
			return Math.abs(greyscale(img,i,j+1)-greyscale(img,i,j-1))+Math.abs(greyscale(img,i-1,j)-greyscale(img,i,j+1));
		}
		else{
			System.out.println("error in c()");
			return 0;
		}
	}
	
	public static double[][] dynamicForwardEnergyMap(double[][][] costs){
		int h = costs.length;
		int w=costs[0].length;
		double[][] M = new double[h][w];
		int i=0;
		int j=0;
		//Upper edges has no forward cost
				M[0][0]=0;
				M[0][w-1]=0;
		//handling the first row
				for(j=1;j<w-1;j++){
					M[i][j]=costs[i][j][0];
					
				}
				//handling the first column (no Cost-left)
				j=0;
				for(i=1;i<h;i++){
					M[i][j]=Math.min(M[i-1][j]+costs[i][j][1], M[i-1][j+1]+costs[i][j][2]);

				}
				//handling the last column (no Cost-right)
				j=w-1;
				for(i=1;i<h;i++){
					M[i][j]=Math.min(M[i-1][j]+costs[i][j][1],M[i-1][j-1]+costs[i][j][0]);
				}
				// handling the rest of the matrix
				for(i=1;i<h;i++){
					for(j=1;j<w-1;j++){
						M[i][j] = Math.min(Math.min(M[i-1][j-1]+costs[i][j][0], M[i-1][j]+costs[i][j][1]),M[i-1][j+1]+costs[i][j][2]);
					
					}
					}
				return M;
	}
	public static BufferedImage readImage(String dest){
		File f=new File(dest);
		BufferedImage img=null;
	double[] times = new double[10];
		try {
			 img = ImageIO.read(f);
		} catch (IOException e) {
			
			e.printStackTrace();
		}
		return img;
	}

	public static void exportEmptyPicture(){
		try{
		    BufferedImage image = null;
			image = new BufferedImage(10, 10, BufferedImage.TYPE_INT_ARGB);
			setPixel( image,0, 0, 30,0);
			setPixel( image,5, 3, 255,99);
			setPixel( image,9, 7, 100,50);
		      File f = new File("images\\Empty.jpg");  //output file path
		      ImageIO.write(image, "jpg", f);
		      System.out.println("Writing complete.");
		    }catch(IOException e){
		      System.out.println("Error: "+e);
		    }
	}


	public static double[][] imageToGreyScaleMatrix(BufferedImage img){
		
		double[][] res = new double[img.getHeight()][img.getWidth()];
		for(int i=0;i<img.getHeight();i++){
			for(int j=0;j<img.getWidth();j++){
				res[i][j] = greyscale(img,j,i);
			}
		}
		return res;
	}
	public static void printMatrix(double[][] m){
		int h = m.length;
		int w = m[0].length;
		for(int i=0;i<h;i++){
			for(int j=0;j<w;j++){
				System.out.print(m[i][j]+"\t");
			}
			System.out.println();
		}
	}
	public static void exportMatrix(double[][] m,String dest){
		PrintWriter out= null;
		try{
		 out = new PrintWriter(dest);}
		catch(Exception e){
			
		}
		
		int h = m.length;
		int w = m[0].length;
		for(int i=0;i<h;i++){
			for(int j=0;j<w;j++){
				out.print((m[i][j]+"\t"));
				
			}
			out.println();
		}
	}
	public static void ExportMatrixAsImage(double [][] matrix,BufferedImage originalImage,String outputDest){
		/**author: Roee
		 * originalImage is need only for the transparency or A value of the pixel
		 * we assume the matrix represent the TRANSPOSE picture - that explains the size match check**/
		int rows= matrix.length;
		int columns = matrix[0].length;
		if(originalImage.getHeight()!= columns || originalImage.getWidth()!= rows){
			System.out.println("Error: unmatched sizes!");
			return;
		}
		BufferedImage image = null;
	    File f = new File(outputDest);
	    image = convertMatrixToImage(matrix,originalImage);
	    try{
		      ImageIO.write(image, "jpg", f);
		    }catch(IOException e){
		      System.out.println("Error: "+e);
		    }
		
	}
	public static BufferedImage convertMatrixToImage(double [][] matrix,BufferedImage originalImage){
		/**author: Roee**/
		int rows= matrix.length;
		int columns = matrix[0].length;
		BufferedImage image = null;
		image = new BufferedImage(originalImage.getWidth(), originalImage.getHeight()
				, BufferedImage.TYPE_INT_ARGB);
		for(int i=0;i<rows;i++){
			for(int j=0;j<columns;j++){
				int pixelVal = (int)Math.round(matrix[i][j]);
				int originalPixelA = (int)Math.round(extractRGB(originalImage,i,j)[3]);
				setPixel(image,i,j,pixelVal,originalPixelA);
				
			}
		}
		return image;
	}
	public static void setPixel(BufferedImage img,int i, int j, int RGB,int A){
		/**author: Roee**/
		int p;
		p = (A<<24) | (RGB<<16) | (RGB<<8) | RGB;
        img.setRGB(i, j, p);	
	}

	public static double[] extractRGB(BufferedImage img,int x, int y){
		 int  clr1   = img.getRGB(x, y);
		 int a1 = (clr1>>24)&0xff;
		 int  r1   = (clr1 & 0x00ff0000) >> 16;
	     int  g1 = (clr1 & 0x0000ff00) >> 8;
	     int  b1  =  clr1 & 0x000000ff;
	     double[] rgb={r1,g1,b1,a1};
	     return rgb;
	     
	}
	public static void setPixel(BufferedImage img,int i, int j, int RGB,int A,int print){
		System.out.println("execute setPixel\t i="+i+"\tj="+j);
		setPixel(img,i,j,RGB,A);
	}

	public static double[] extractRGB(BufferedImage img,int x, int y,int print){
		System.out.println("execute extractRGB\t x="+x+"\ty="+y);
		return extractRGB(img,x, y);
	     
	}

	public static double[][] addEntropy(double[][] Emap,BufferedImage OrigImg,int memoizationInd){
		/**author: Roee
		 * part 1 - 2
		 * memoizationInd = 0 for no memoization, 1 for P-Values memo,
		 *  2 for greyscale memom, 3 for both memo**/
		HashMap<Integer,HashMap<Integer,Double>> pValuesDictionary = new HashMap<>();
		HashMap<Integer,HashMap<Integer,Double>> greyscaleDictionary = new HashMap<>();
		int rows = Emap.length;
		int columns = Emap[0].length;
		double[][] res = new double[rows][columns];
		for(int i=0;i<rows;i++){
			for(int j=0;j<columns;j++){
				res[i][j]=Emap[i][j]+pixelEntropy(OrigImg,i,j,
						pValuesDictionary,greyscaleDictionary,memoizationInd);
			}
		}
		
		
		return res;
	}
	
	public static double[][] returnEntropy(BufferedImage OrigImg,int memoizationInd){
		/**author: Roee
		 * part 1 - 2
		 * memoizationInd = 0 for no memoization, 1 for P-Values memo,
		 *  2 for greyscale memom, 3 for both memo**/
		HashMap<Integer,HashMap<Integer,Double>> pValuesDictionary = new HashMap<>();
		HashMap<Integer,HashMap<Integer,Double>> greyscaleDictionary = new HashMap<>();
		int rows = OrigImg.getHeight();
		int columns = OrigImg.getWidth();
		double[][] res = new double[rows][columns];
		for(int i=0;i<rows;i++){
			for(int j=0;j<columns;j++){
				res[i][j]=pixelEntropy(OrigImg,i,j,
						pValuesDictionary,greyscaleDictionary,memoizationInd);
			}
		}
		
		
		return res;
	}
	public static double pixelEntropy(BufferedImage OrigImg,int i,int j,
			HashMap<Integer,HashMap<Integer,Double>> pValuesDictionary,
			HashMap<Integer,HashMap<Integer,Double>> greyscaleDictionary,int memoizationInd){
		int rows = OrigImg.getHeight();
		int columns = OrigImg.getWidth();

		int left = Math.max(0, i-4);
		int right = Math.min(columns-1, i+4);
		int up = Math.max(0, j-4);
		int down = Math.min(rows-1, j+4);
		double sum=0;
		double p;
		double plogp;
		int neighbors=0;
		
		for(int m=left;m<=right;m++){
			for(int n=up;n<=down;n++){
				neighbors++;
				if(memoizationInd==1||memoizationInd==3){
					if(dictionaryContainsKey(pValuesDictionary,m,n)){
						plogp= dictionaryGet(pValuesDictionary,m,n);
						Pcounter++;
					}
					else{
						p=p(OrigImg,m,n,greyscaleDictionary,memoizationInd);
						plogp=p*Math.log(p);
						dictionaryPut(pValuesDictionary,m,n,plogp);
						}
					}
				else{
					p=p(OrigImg,m,n,greyscaleDictionary,memoizationInd);
					plogp=p*Math.log(p);
					
					}
					sum+=plogp;
					
				
				}	
		}
		return -sum/neighbors;
	}
	/*
	public static double pixelEntropySlidingWindow(BufferedImage OrigImg,int i,int j,
			HashMap<Integer,HashMap<Integer,Double>> pValuesDictionary,
			HashMap<Integer,HashMap<Integer,Double>> greyscaleDictionary,int memoizationInd,double prevEntropy){
		int rows = OrigImg.getHeight();
		int columns = OrigImg.getWidth();

		int left = Math.max(0, i-4);
		int right = Math.min(columns-1, i+4);
		int up = Math.max(0, j-4);
		int down = Math.min(rows-1, j+4);
		double leftColumnOfPrev = calculateColumnEntropy(i,j,-5,)
		double rightColumnOfCurrent;
		
		
		
		
	}
	public static double calculateColumnEntropy(int i;int j;int columnShift;){
		
	}*/

	public static double p(BufferedImage OrigImg,int m,int n,
			HashMap<Integer,
			HashMap<Integer,Double>> greyscaleDictionary,int memoizationInd){
		int rows = OrigImg.getHeight();
		int columns = OrigImg.getWidth();
		int left = Math.max(0, m-4);
		int right = Math.min(columns-1, m+4);
		int up = Math.max(0, n-4);
		int down = Math.min(rows-1, n+4);
		double greyscale;
		double sum=0;
		int neighbors=0;
		for(int k=left;k<=right;k++){
			for(int l=up;l<=down;l++){
				neighbors++;
				if(memoizationInd==2||memoizationInd==3){
					if(dictionaryContainsKey(greyscaleDictionary,k,l)){
						greyscale= dictionaryGet(greyscaleDictionary,k,l);
						Gcounter++;
					}
					else{
						greyscale=greyscale(OrigImg,k,l);
						dictionaryPut(greyscaleDictionary,k,l,greyscale);
				}
				}
				else{
					greyscale=greyscale(OrigImg,k,l);
				}
				
				sum+=greyscale;
		}}
		
		return greyscale(OrigImg,m,n)/(sum/neighbors);
	}
	public static double greyscale(BufferedImage OrigImg,int m,int n){
		double[] RGBA = extractRGB(OrigImg,m, n);
		double avg = (RGBA[0]+RGBA[1]+RGBA[2])/3;
		return avg;
	}
	public static void dictionaryPut(HashMap<Integer,HashMap<Integer,Double>> d,
			 int i, int j,double val){
		 if(d.get(i)==null){
			 d.put(i, new HashMap<Integer,Double>());
		 }
		 d.get(i).put(j, val);
	 }
	 public static double dictionaryGet(HashMap<Integer,HashMap<Integer,Double>> d,
			 int i, int j){
		 return d.get(i).get(j);
	 }
	 public static boolean dictionaryContainsKey(HashMap<Integer,HashMap<Integer,Double>> d,
			 int i, int j){
		 if(d.containsKey(i)){
			 if(d.get(i).containsKey(j)){
				 return true;
			 }
		 }
		 return false;
		 
		 

		 
	 }
	public static double[][] dynamicMap( double[][] Emap,int direction){
		double[][] dyMap=new double[Emap.length][Emap[0].length];
		if(direction==0) {// for horizontal
			Emap=transposeMatrix(Emap);
			dyMap=transposeMatrix(dyMap);
			
			
				}
		
		for(int i=0;i<Emap.length;i++) {
			for (int j=0; j<Emap[0].length; j++) {
				if (i==0) {
					dyMap[i][j]=Emap[i][j];
					
					
				}
				else if(j!=0 && j!=Emap[0].length-1) {
					double minPath=Math.min(dyMap[i-1][j-1] , dyMap[i-1][j]) ;
					minPath=Math.min(minPath   ,  dyMap[i-1][j+1]);
					dyMap[i][j]=Emap[i][j]+minPath;
					
					
				}
				else if(j==0) {
					dyMap[i][j]=Emap[i][j]+Math.min( dyMap[i-1][j]  ,dyMap[i-1][j+1]) ;
				}
				
				else if(j==Emap[0].length-1) {
					dyMap[i][j]=Emap[i][j]+Math.min( dyMap[i-1][j-1]  ,dyMap[i-1][j]) ;
				}
			}
			/*if(direction==1) {
				Emap=transposeMatrix(Emap);// transpose image again in case it was horizontal
				dyMap=transposeMatrix(dyMap);// transpose image again in case it was vertical- might be REDUNDANT!
			
					}*/
		}
		if(direction==0) {
			Emap=transposeMatrix(Emap);// transpose image again in case it was horizontal
			dyMap=transposeMatrix(dyMap);// transpose image again in case it was vertical- might be REDUNDANT!
		
				}
		
		
		
		return dyMap;
	}
	 public static double[][] transposeMatrix(double [][] m){ // THIS WAS ADDED 
	        double[][] temp = new double[m[0].length][m.length];
	        for (int i = 0; i < m.length; i++)
	            for (int j = 0; j < m[0].length; j++)
	                temp[j][i] = m[i][j];
	        return temp;
	    }
	public static int[] chooseSeam(double[][] Dmap, int direction){
		/* part 1 - 4
		 * direction = 0 for horizontal, 1 for verticle
		 * 
		 * */
		int height=Dmap[0].length;
		int width=Dmap.length;
		int seam[] = null;
		if(direction ==0) {// horizontal seam - the seam includes for each column y its only x which is in seam
			 seam=new int[height];
			double first_min=Dmap[0][height - 1];
			int index=0;
			
			for(int x=0;x<width;x++) {
				if(Dmap[x][height - 1]<first_min) {
					
					first_min=Dmap[x][height - 1];
					index=x;
				}
				
				
			}
			seam[height-1]=index;
		
			int x=index;
			for(int y=height-2;y>=0;y--) {
				if(x!=0 && x!=width-1 ) {
				double min=Math.min(Dmap[x][y],Dmap[x+1][y]);
				min=Math.min(min, Dmap[x-1][y]);
				if (Dmap[x][y]== min)
					seam[y]=x;
				if (Dmap[x-1][y]== min) {
					seam[y]=x-1;
					x=x-1;
				}
				if (Dmap[x+1][y]== min) {
					seam[y]=x+1;
					x=x+1;}
				
			}
				else if(x==0) {
					double min=Math.min(Dmap[x][y],Dmap[x+1][y]);
					if (Dmap[x+1][y]== min) {
						seam[y]=x+1;
						x=x+1;}
					if (Dmap[x][y]== min)
						seam[y]=x;
					
				}
				else if(x==width-1) {
					double min=Math.min(Dmap[x][y],Dmap[x-1][y]);
					if (Dmap[x-1][y]== min) {
						seam[y]=x-1;
						x=x-1;}
					if (Dmap[x][y]== min)
						seam[y]=x;
					
				}
			}
			
		}
		if(direction ==1) {// vertical seam- for each row in map includes its only y in seam
			seam=new int[width];
			double first_min=Dmap[width-1][0];// bottom row
			int index=0;
			
			for(int y=0;y<height;y++) {
				if(Dmap[width-1][y]<first_min) {
					
					first_min=Dmap[width-1][y];
					index=y;
				}
				
				
			}
			seam[width-1]=index;
			int y=index;
			for(int x=width-2;x>=0;x--) {
				if(y!=0 && y!=height-1) {
					double min=Math.min(Dmap[x][y+1],Dmap[x][y]);
					min=Math.min(min, Dmap[x][y-1]);
					if (Dmap[x][y]== min)
						seam[x]=y;
					if (Dmap[x][y-1]== min) {
						seam[x]=y-1;
						y=y-1;
					}
					if (Dmap[x][y+1]== min) {
						seam[x]=y+1;
						y=y+1;}
				}
				else if(y==0) {
					double min=Math.min(Dmap[x][y],Dmap[x][y+1]);
					if (Dmap[x][y+1]== min) {
						seam[x]=y+1;
						y=y+1;}
					if (Dmap[x][y]== min)
						seam[x]=y;
					
				}
				else if(y==height-1) {
					double min=Math.min(Dmap[x][y],Dmap[x][y-1]);
					if (Dmap[x][y-1]== min) {
						seam[x]=y-1;
						y=y-1;}
					if (Dmap[x][y]== min)
						seam[x]=y;
					
				}
				
			}
			
			
		}
		return seam;
	}
	public static BufferedImage deleteSeam(BufferedImage img, int direction, int[] seamVector){
		/* part 1 - 4
		 * direction = 0 for vertical, 1 for horizontal
		 * 
		 * */
		int width = img.getWidth();
		int height = img.getHeight();
		BufferedImage bufferedImage = null;

		if(direction==0) {// vertical seam
			 bufferedImage = new BufferedImage(width-1, height, BufferedImage.TYPE_INT_RGB);
			boolean shift;
			for(int y=0; y<height; y++) {
				 shift=false;
				for(int x=0;x<width;x++) {
					if (seamVector[y]==x) {
						shift=true;
					}
					else {
						if (shift==true) {
							bufferedImage.setRGB(x-1, y,img.getRGB(x,y));
						}
						else {
						bufferedImage.setRGB(x, y,img.getRGB(x,y));
						}
					}
					
				}
			}
			
			
		}
		else if(direction==1) {// cut horizontal
			 bufferedImage = new BufferedImage(width, height-1, BufferedImage.TYPE_INT_RGB);
			boolean shift;
			for(int x=0; x<width; x++) {
				 shift=false;
				for(int y=0;y<height;y++) {
					if (seamVector[x]==y) {
						shift=true;
						
					}
					else {
						if (shift==true) {
							bufferedImage.setRGB(x, y-1,img.getRGB(x,y));
						}
						else {
						bufferedImage.setRGB(x, y,img.getRGB(x,y));
						}
					}
					
				}
			}
			
		}
		return bufferedImage;
	}

	public static double[][] createEnegryMap(BufferedImage img,int width,int height){
		double energy_map[][] = new double[width][ height];
		
		for (int x = 0; x < width; x++) {
		    for (int y = 0; y < height; y++) {
		    	energy_map[x][y]=energy(x,y,img);
		    } 
		}
		return energy_map;
		
	}

	public static double energy(int x,int y,BufferedImage img) {
		double energy=0;
		int width = img.getWidth();
		int height = img.getHeight();
		
		 if(x>0 && y<height-1&& x<width-1 &&y>0) {
		   energy=(val(x,y,x-1,y-1,img)+val(x,y,x-1,y,img)+val(x,y,x-1,y+1,img)
		            +val(x,y,x,y-1,img)+val(x,y,x,y+1,img)+
				   val(x,y,x+1,y-1,img)+val(x,y,x+1,y,img)+val(x,y,x+1,y+1,img) ) /8;
				   
		}
		 else if(x==0 &&y==0 ) {
			 energy=(val(x,y,x,y+1,img)+val(x,y,x+1,y,img)+val(x,y,x+1,y+1,img) ) /3; 
		 }
		 
		 else if( x==width-1 && y==0  ) {
			 energy=(val(x,y,x-1,y,img)+val(x,y,x,y+1,img)+val(x,y,x-1,y+1,img) ) /3; 
		 }
		 
		 else if( x==0 && y==height-1) {
			 energy=(val(x,y,x+1,y,img)+val(x,y,x+1,y-1,img)+val(x,y,x,y-1,img) ) /3; 
		 }
		 
		 else if( x==width-1  && y==height-1) {
			 energy=(val(x,y,x,y-1,img)+val(x,y,x-1,y,img)+val(x,y,x-1,y-1,img) ) /3; 
		 }
		 
		 else if(x== width-1) {
			  energy=(val(x,y,x-1,y-1,img)+val(x,y,x-1,y,img)+val(x,y,x-1,y+1,img)+val(x,y,x,y-1,img)+val(x,y,x,y+1,img)) /5;
		 }
		 
		 else if(y==height-1){
			 energy=(val(x,y,x+1,y,img)+val(x,y,x+1,y-1,img)+val(x,y,x,y-1,img)+val(x,y,x-1,y,img)+val(x,y,x-1,y-1,img)) /5;
		 }
		 
		 else if(y==0){
			 energy=(val(x,y,x-1,y,img)+val(x,y,x-1,y+1,img)+val(x,y,x,y+1,img)+val(x,y,x+1,y+1,img)+val(x,y,x+1,y,img)) /5;
		 }
		 
		 else if(x==0){
			 energy=(val(x,y,x,y+1,img)+val(x,y,x+1,y+1,img)+val(x,y,x+1,y,img)+val(x,y,x+1,y-1,img)+val(x,y,x,y-1,img)) /5;
		 }
		 return energy;
		 
		
	}

	private static double val(int x, int y, int i, int j,BufferedImage img) {
	
        double r1=extractRGB(img,x,y)[0];
        double g1=extractRGB(img,x,y)[1];
        double b1=extractRGB(img,x,y)[2];
        
        double r2=extractRGB(img,i,j)[0];
        double g2=extractRGB(img,i,j)[1];
        double b2=extractRGB(img,i,j)[2];
        
		return (Math.abs(r1-r2)+ Math.abs(g1-g2) +Math.abs(b1-b2))/3;
	}


}

