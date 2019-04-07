
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
		//exportDynamicMaps("strawberry",0.5);
		//entropyCorrectnessTest("surfer");
		//chooseAndDeleteSeamTest("baloon",250,0,0);
		//exportEmptyPicture();
		forwardEnergyTest("cat",100 ,0.0);
		/*forwardEnergyTest("cat",100 ,0.5);
		forwardEnergyTest("cat",100 ,0.8);
		forwardEnergyTest("cat",100 ,1.0);*/
		

	}
	
	public static void chooseAndDeleteSeamTest(String fileName,int seams,int direction,double weight){
		BufferedImage img = readImage("images\\"+fileName+".jpg");
		double[] times = new double[10];
		int absEntropy=1;
		int scalar = 50;
		times[0] = System.nanoTime()/1000000000;
		for(int i=0;i<seams;i++){
			System.out.println(i);
			double[][] Emap = createEnegryMap(img, img.getWidth(),img.getHeight());
			double[][] entropy = returnEntropySlidingWindow(img,1,0);
			if(absEntropy==1){
				matrixAbsoloute(entropy);}
			multiplyMatrixByScalar(entropy,scalar);
			Emap = weightedAverageOfMatrices(entropy,Emap,weight);
			double[][] Dmap = dynamicMap(Emap,direction);
			int[] seam =  chooseSeam(Dmap,direction);
			img= deleteSeam( img, direction, seam);
			}
		times[1] = System.nanoTime()/1000000000;
		System.out.println("total time: "+(times[1]-times[0]));
		File f = new File("images\\"+fileName+"\\w="+weight+"removed "+seams+"seams"+" direction="+direction+"absEntropy="+absEntropy+"scalarEntropy="+scalar+".jpg");
		try{
	    ImageIO.write(img, "jpg", f);}
		catch(IOException e){
		      System.out.println("Error: "+e);
		    }
		
		
	}
	
	public static void entropyCorrectnessTest(String fileName){
		BufferedImage image = null;
		image = readImage("images\\"+fileName+".jpg");
		int width = image.getWidth();
		int height = image.getHeight();
		int scalar=100;
		double w=0.5;
		double eMap[][] = createEnegryMap(image, width,height);
		exportMatrixToTextFile(eMap,"eMap.txt");
		double[][] Entmap = returnEntropySlidingWindow(image,1,1);
		//matrixAbsoloute(Entmap);
		multiplyMatrixByScalar(Entmap,scalar);
		double[] energyRange = checkValuesRange(eMap);
		double[] EntmapRange = checkValuesRange(Entmap);
		System.out.println("energy:"+energyRange[0]+","+energyRange[1]+","+energyRange[2]);
		System.out.println("Entmap:"+EntmapRange[0]+","+EntmapRange[1]+","+EntmapRange[2]);
		
		double[][] weighted = weightedAverageOfMatrices(Entmap,eMap,w);
		ExportMatrixAsImage(eMap,image,"images\\"+fileName+"\\eMap.jpg");
		ExportMatrixAsImage(Entmap,image,"images\\"+fileName+"\\entropy normalizedX+"+scalar+".jpg");
		ExportMatrixAsImage(weighted,image,"images\\"+fileName+"\\weighted w="+w+".jpg");
		chooseAndDeleteSeamTest("surfer",50,0,w);
		//exportMatrixToTextFile(Entmap,"EntmapX"+scalar+".txt");
	}
	public static void multiplyMatrixByScalar(double[][] m,int s){
		int r = m.length;
		int c = m[0].length;
		
		for(int i=0;i<r;i++){
			for(int j=0;j<c;j++){
				m[i][j] = m[i][j]*s;
			}
		}
		
	}
	public static void matrixAbsoloute(double[][] m){
		int r = m.length;
		int c = m[0].length;
		
		for(int i=0;i<r;i++){
			for(int j=0;j<c;j++){
				m[i][j] = Math.abs(m[i][j]);
			}
		}
		
	}
public static double[][] transpose(double[][] m){
	int r = m.length;
	int c = m[0].length;
	double[][] res = new double[c][r];
	for(int i=0;i<r;i++){
		for(int j=0;j<c;j++){
			res[j][j] = m[i][j];
		}
	}
	return res;
		
	}
	public static double[] checkValuesRange(double[][] m){
		/**for testing only - to check the wanted weight between engery map and entropy
		 * res[0] = min, res[1] = max, res[2] = avg**/
		int r = m.length;
		int c = m[0].length;
		double min=1.0/0.0;//for positive infinity
		double max=-1.0/0.0;//for negative infinity
		double sum = 0;
		double v;
		for(int i=0;i<r;i++){
			for(int j=0;j<c;j++){
				v = m[i][j];
				sum+=v;
				if(v>max){
					max = v;
				}
				if(v<min){
					min = v;
				}
			}
		}
			double[] res=new double[3];
			res[0]=min;
			res[1]=max;
			res[2]=sum/((r+1)*(c+1));
			return res;
		
		
	}

	public static void exportEmptyPicture(){
		int size=400;
		int halfsize = Math.round(size/2);
		try{
		    BufferedImage image = null;
			image = new BufferedImage(size, size, BufferedImage.TYPE_INT_ARGB);
			for(int i=0;i<halfsize;i++){
				for(int j=0;j<size;j++){
					setPixel( image,i, j, -10,255);
				}
			}
		      File f = new File("images\\Empty.jpg");  //output file path
		      ImageIO.write(image, "jpg", f);
		      System.out.println("Writing complete.");
		    }catch(IOException e){
		      System.out.println("Error: "+e);
		    }
	}

	public static int[] areIdenticMatrices(double[][] m1,double[][]m2){
		int rows=m1.length;
		int columns = m1[0].length;
		if(m2.length!=rows||m2[0].length!=columns){
			System.out.println("Error unmatch sizes");
			return null;
		}
		for(int i=0;i<rows;i++){
			for(int j=0;j<columns;j++){
				double val1 = m1[i][j];
				double val2 = m2[i][j];
				if(Math.abs(val1-val2)>0.000000001){
					int[] res = {i,j};
					return res;
				}
			}
		}
		int[] res = {-1,-1};
		return res;
	}
	public static void entropyTimeTest(String fileName){
		int[] memoizationValues= {1,3};
		 File f=new File("images\\"+fileName+".jpg");
			BufferedImage img=null;
		double[] times = new double[10];
			try {
				 img = ImageIO.read(f);
			} catch (IOException e) {
				
				e.printStackTrace();
			}
			
			int width = img.getWidth();
			int height = img.getHeight();
			System.out.println("image Height:"+img.getHeight()+"image Width:"+img.getWidth());
		
			times[0] = System.nanoTime()/1000000000;
			double energy_map[][] = createEnegryMap(img, width,height);
			System.out.println("Emap Height:"+energy_map.length+"Emap Width:"+energy_map[0].length);
		
			//System.out.println("finish creating energy map");
			times[1] = System.nanoTime()/1000000000;
			System.out.println("energy time:"+(times[1]-times[0]));
			ExportMatrixAsImage(energy_map,img,"images\\"+fileName+"\\"+fileName+"energyMap.jpg");
				for(int p=0;p<2;p++){
					times[2] = System.nanoTime()/1000000000;
					double[][] energy_map_after_entropy =addEntropySlidingWindow(energy_map,img,memoizationValues[p],0);
					times[3] = System.nanoTime()/1000000000;
					System.out.println("entorpy time:"+(times[3]-times[2]));
					ExportMatrixAsImage(energy_map_after_entropy,img,"images\\"+fileName+"\\"+fileName+" enrgyWithEntropy.jpg");
					System.out.println("memoization Ind:"+memoizationValues[p]);
					System.out.println("Pcounter:"+Pcounter);
					System.out.println("Gcounter:"+Gcounter+"\n");
					Pcounter=Gcounter=0;
					}
				
	}
	public static void forwardEnergyTest(String fileName,int seams,double weight){
		
		BufferedImage img = readImage("images\\"+fileName+".jpg");
		for(int i=0;i<seams;i++){
			System.out.println("i= "+i+"\tw= "+weight);
			double[][][] costs = forwardEnergyCost(img);
			double[][] dynamicForwardMap = dynamicForwardEnergyMap(costs);
			double[][] Emap = createEnegryMap(img, img.getWidth(),img.getHeight());
			double[][] dynamicEnergyMap = dynamicMap(Emap,0);
			double [][]Dmap = weightedAverageOfMatrices(dynamicForwardMap,dynamicEnergyMap,weight);
			int[] seam =  chooseSeam(Dmap,0);
			img= deleteSeam( img, 0, seam);
		
		}
		
			
		
		File f = new File("images\\Forward: "+fileName+" w="+weight+"removed "+seams+"seams"+".jpg");
		try{
	    ImageIO.write(img, "jpg", f);}
		catch(IOException e){
		      System.out.println("Error: "+e);}
	}
	public static void exportDynamicMaps(String fileName,double w){
		BufferedImage img = readImage("images\\"+fileName+".jpg");
		double[][][] costs = forwardEnergyCost(img);
		double[][] dynamicForwardMap = dynamicForwardEnergyMap(costs);
		double[][] Emap = createEnegryMap(img, img.getWidth(),img.getHeight());
		double[][] dynamicEnergyMap = dynamicMap(Emap,0);
		double [][]weightedMap = weightedAverageOfMatrices(dynamicForwardMap,dynamicEnergyMap,w);
		ExportMatrixAsImage(averageOfCosts(costs),img,"images\\"+fileName+"\\"+" costs.jpg");
		//ExportMatrixAsImage(dynamicEnergyMap,img,"images\\"+fileName+"\\"+ "dynamicEnergyMap.jpg");
		//ExportMatrixAsImage(Emap,img,"images\\"+fileName+"\\"+ "EnergyMap.jpg");
		ExportMatrixAsImage(dynamicForwardMap,img,"images\\"+fileName+"\\"+ "dynamicForwardMap.jpg");
		//ExportMatrixAsImage(weightedMap,img,"images\\"+fileName+"\\"+ "weightedMap"+w+".jpg");
		
		
		exportMatrixToTextFile(dynamicForwardMap, fileName+"dynamicForward.txt");
		exportMatrixToTextFile(averageOfCosts(costs), fileName+"averageCosts.txt");
		exportMatrixToTextFile(dynamicEnergyMap, fileName+"dynamicEnergyMap.txt");
	}
	public static double[][] averageOfCosts(double[][][] costs){
		int n=costs.length;
		int m = costs[0].length;
		double[][] res = new double[n][m];
		for(int i=0;i<n;i++){
			for(int j=0;j<m;j++){
				res[i][j] = (costs[i][j][0]+costs[i][j][1]+costs[i][j][2])/3;
			}
		}
		return res;
		
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
			res[i][j][0]=Math.abs(oldGreyScaleForForwarding(img,i,j-1)-oldGreyScaleForForwarding(img,i,j-1));

			
		}

		//handling the first column (no Cost-left)
		j=0;
		for(i=1;i<h;i++){

			res[i][j][1]=0;
			res[i][j][2]=Math.abs(oldGreyScaleForForwarding(img,i-1,j)-oldGreyScaleForForwarding(img,i,j+1));
		}
		//handling the last column (no Cost-right)
		j=w-1;
		for(i=1;i<h;i++){
			res[i][j][1]=0;
			res[i][j][0]=Math.abs(oldGreyScaleForForwarding(img,i-1,j)-oldGreyScaleForForwarding(img,i,j-1));
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
			return Math.abs(oldGreyScaleForForwarding(img,i,j+1)-oldGreyScaleForForwarding(img,i,j-1))+Math.abs(oldGreyScaleForForwarding(img,i-1,j)-oldGreyScaleForForwarding(img,i,j-1));
		}
		if(direction=='U'){
			return Math.abs(oldGreyScaleForForwarding(img,i,j+1)-oldGreyScaleForForwarding(img,i,j-1));
		}
		if(direction=='R'){
			return Math.abs(oldGreyScaleForForwarding(img,i,j+1)-oldGreyScaleForForwarding(img,i,j-1))+Math.abs(oldGreyScaleForForwarding(img,i-1,j)-oldGreyScaleForForwarding(img,i,j+1));
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
	public static void exportMatrixToTextFile(double[][] m,String dest){
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
		int pixelVal = 0;
		BufferedImage image = null;
		image = new BufferedImage(originalImage.getWidth(), originalImage.getHeight()
				, BufferedImage.TYPE_INT_ARGB);
		for(int i=0;i<rows;i++){
			for(int j=0;j<columns;j++){
				pixelVal = (int)Math.round(matrix[i][j]);
				int originalPixelA = (int)Math.round(extractRGB(originalImage,i,j)[3]);

				setPixel(image,i,j,pixelVal,255);//I've change the originalPixelA to 255 to see the full map
				
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
	public static double[][] returnEntropy(BufferedImage OrigImg,int memoizationInd,int normalizationInd){
		/**author: Roee
		 * return the transpose of the Entropy - because we want it to suit the transposed Energy-map
		 * memoizationInd = 0 for no memoization, 1 for P-Values memo,
		 *  2 for greyscale memom, 3 for both memo**/
		//rows and columns are revered because we want to return the transpose of the Entropy
		int rows = OrigImg.getWidth();
		int columns = OrigImg.getHeight();
		double[][] res = new double[rows][columns];
		HashMap<Integer,HashMap<Integer,Double>> pValuesDictionary = new HashMap<>();
		HashMap<Integer,HashMap<Integer,Double>> greyscaleDictionary = new HashMap<>();
		for(int i=0;i<rows;i++){
			for(int j=0;j<columns;j++){
				res[i][j]=pixelEntropy(OrigImg,j,i,//i,j are reversed the order because the Emap is transposed
						pValuesDictionary,greyscaleDictionary,memoizationInd);
			}
		}
		if(normalizationInd==1){
			res = normalize(res);
		}
		
		return res;
	}

	public static double[][] addEntropy(double[][] Emap,BufferedImage OrigImg,int memoizationInd,int normalizationInd){
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
				res[i][j]=Emap[i][j]+pixelEntropy(OrigImg,j,i,//i,j are reversed the order because the Emap is transposed
						pValuesDictionary,greyscaleDictionary,memoizationInd);
			}
		}
		if(normalizationInd==1){
			res =normalize(res);
		}
		return res;
	}
	
	public static double[][] addEntropySlidingWindow(double[][] Emap,BufferedImage OrigImg,int memoizationInd,int normalizationInd){
		/**author: Roee
		 * part 1 - 2
		 * memoizationInd = 0 for no memoization, 1 for P-Values memo,
		 *  2 for greyscale memom, 3 for both memo**/
		HashMap<Integer,HashMap<Integer,Double>> pValuesDictionary = new HashMap<>();
		HashMap<Integer,HashMap<Integer,Double>> greyscaleDictionary = new HashMap<>();
		int rows = Emap.length;
		int columns = Emap[0].length;
		double pixelEntropy;
		int j;
		double[][] res = new double[rows][columns];
		for(int i=0;i<rows;i++){
			j=0;
			pixelEntropy=pixelEntropy(OrigImg,j,i,/*i,j are reversed the order because the Emap is transposed*/
					pValuesDictionary,greyscaleDictionary,memoizationInd);
			res[i][j]=Emap[i][j]+pixelEntropy;
			for(j=1;j<columns;j++){
				pixelEntropy=pixelEntropySlidingWindow(OrigImg,j,i,/*i,j are reversed the order because the Emap is transposed*/
						pValuesDictionary,greyscaleDictionary,memoizationInd,pixelEntropy);
				res[i][j]=Emap[i][j]+pixelEntropy;
			}
		}
		if(normalizationInd==1){
			res = normalize(res);
		}
		return res;
	}
	public static double[][] returnEntropySlidingWindow(BufferedImage OrigImg,int memoizationInd,int normalizationInd){
		/**author: Roee
		 * part 1 - 2
		 * memoizationInd = 0 for no memoization, 1 for P-Values memo,
		 *  2 for greyscale memom, 3 for both memo**/
		HashMap<Integer,HashMap<Integer,Double>> pValuesDictionary = new HashMap<>();
		HashMap<Integer,HashMap<Integer,Double>> greyscaleDictionary = new HashMap<>();
		int rows = OrigImg.getWidth();
		int columns = OrigImg.getHeight();
		double pixelEntropy;
		int j;
		double[][] res = new double[rows][columns];
		for(int i=0;i<rows;i++){
			j=0;
			pixelEntropy=pixelEntropy(OrigImg,j,i,/*i,j are reversed the order because the Emap is transposed*/
					pValuesDictionary,greyscaleDictionary,memoizationInd);
			res[i][j]=pixelEntropy;
			for(j=1;j<columns;j++){
				pixelEntropy=pixelEntropySlidingWindow(OrigImg,j,i,/*i,j are reversed the order because the Emap is transposed*/
						pValuesDictionary,greyscaleDictionary,memoizationInd,pixelEntropy);
				res[i][j]=pixelEntropy;
			}
		}
		if(normalizationInd==1){
			res = normalize(res);
		}
		return res;
	}
	
	public static double[][] returnPValueMatrix(BufferedImage OrigImg){
		/**for testing the entropy only
		 * returns the TRANSPOSED P-values matrix of the image**/
		//rows and columns are revered because we want to return the transpose of the Entropy
		int rows = OrigImg.getWidth();
		int columns = OrigImg.getHeight();
		double[][] res = new double[rows][columns];
		for(int i=0;i<rows;i++){
			for(int j=0;j<columns;j++){
				res[i][j]=p(OrigImg,j,i,null,0);
			}
		}
		
		
		return res;
	}
	public static double handleNanValues(double x){
		if(Double.isNaN(x)){
			return 0;
		}
		return x;
	}
	public static double pixelEntropy(BufferedImage OrigImg,int i,int j,
			HashMap<Integer,HashMap<Integer,Double>> pValuesDictionary,
			HashMap<Integer,HashMap<Integer,Double>> greyscaleDictionary,int memoizationInd){
		/**returns the entropy of the pixel(i,j) in the picture - i is the row number**/
		//set boundaries part
		int rows = OrigImg.getHeight();
		int columns = OrigImg.getWidth();
		int left = Math.max(0, j-4);
		int right = Math.min(columns-1, j+4);
		int up = Math.max(0, i-4);
		int down = Math.min(rows-1, i+4);
		double sum=0;
		double p;
		double plogp;
		//loop part
		for(int m=up;m<=down;m++){
			for(int n=left;n<=right;n++){

				if(memoizationInd==1||memoizationInd==3){
					if(dictionaryContainsKey(pValuesDictionary,m,n)){
						plogp= dictionaryGet(pValuesDictionary,m,n);
						Pcounter++;
					}
					else{
						p=p(OrigImg,m,n,greyscaleDictionary,memoizationInd);
						plogp=handleNanValues(p*Math.log(p));
						//handleNanValues is needed because p could be 0 in which case plogp=Nan
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
		return -sum;
	}
	public static int countNeighbors(int i,int j,int rows,int columns){
		/**returns the number of neighbors the (i,j) pixel has in
		 * a 9X9 window centered at the pixel in the picture**/
		int left = Math.max(0, j-4);
		int right = Math.min(columns-1, j+4);
		int up = Math.max(0, i-4);
		int down = Math.min(rows-1, i+4);
		return (right-left+1)*(down-up+1);
	}
	public static double[][] normalize(double[][] m){
		int r = m.length;
		int c = m[0].length;
		double[][] res = new double[r][c];
		for(int i=0;i<r;i++){
			for(int j=0;j<c;j++){
				res[i][j] = m[i][j]/countNeighbors(i,j,r,c);
			}
		}
		return res;
	}
	public static double pixelEntropySlidingWindow(BufferedImage OrigImg,int i,int j,
			HashMap<Integer,HashMap<Integer,Double>> pValuesDictionary,
			HashMap<Integer,HashMap<Integer,Double>> greyscaleDictionary,int memoizationInd,double prevEntropy){
		/**@pre: prevEntropy=pixelEntropy(i-1,j)**/
		double sum;
		double upperRowOfPrev = calculateSlidingWindowRowEntropy(OrigImg,i,j,-5,pValuesDictionary,
				greyscaleDictionary, memoizationInd,prevEntropy);//-5 because we need the left column of the PREVIOUS sliding window
		double lowerRowOfCurrent= calculateSlidingWindowRowEntropy(OrigImg,i,j,4,pValuesDictionary,
				greyscaleDictionary, memoizationInd,prevEntropy);
		sum=(-prevEntropy)-upperRowOfPrev+lowerRowOfCurrent;
		//we multiply prevEntropy by (-1) because prevEntropy=-(sum of plogp values)
		return -sum;
		
		
		
	}
	public static double calculateSlidingWindowRowEntropy(BufferedImage OrigImg,int i,int j,int rowShift,
			HashMap<Integer,HashMap<Integer,Double>> pValuesDictionary,
			HashMap<Integer,HashMap<Integer,Double>> greyscaleDictionary,int memoizationInd,double prevEntropy){
		/**calculate the sum of plogp values in the relevant row of the sliding window
		 * 9X9 sliding window**/
		int rows = OrigImg.getHeight();
		int columns = OrigImg.getWidth();
		int m = i+rowShift;
		if(m<0||m>=rows){
			//left and right edges of the picture, there is no such column
			return 0;
		}
		int left = Math.max(0, j-4);
		int right = Math.min(columns-1, j+4);
		double sum =0;
		double plogp;
		double p;
		for(int n=left;n<=right;n++){
			if(memoizationInd==1||memoizationInd==3){
				if(dictionaryContainsKey(pValuesDictionary,m,n)){
					plogp= dictionaryGet(pValuesDictionary,m,n);
					Pcounter++;
				}
				else{
					p=p(OrigImg,m,n,greyscaleDictionary,memoizationInd);
					plogp=handleNanValues(p*Math.log(p));
					//handleNanValues is needed because p could be 0 in which case plogp=Nan
					dictionaryPut(pValuesDictionary,m,n,plogp);
					}
				}
			else{
				p=p(OrigImg,m,n,greyscaleDictionary,memoizationInd);
				plogp=p*Math.log(p);
				
				}
				sum+=plogp;
		}
		return sum;
	}

	public static double p(BufferedImage OrigImg,int m,int n,
			HashMap<Integer,
			HashMap<Integer,Double>> greyscaleDictionary,int memoizationInd){
		/**m=row number, n=columns number (according to the instructions page)**/
		int rows = OrigImg.getHeight();
		int columns = OrigImg.getWidth();
		int left = Math.max(0, n-4);
		int right = Math.min(columns-1, n+4);
		int up = Math.max(0, m-4);
		int down = Math.min(rows-1, m+4);
		double greyscale;
		double sum=0;
		int neighbors=0;
		for(int k=up;k<=down;k++){
			for(int l=left;l<=right;l++){
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
		
		return handleNanValues(greyscale(OrigImg,m,n)/(sum/neighbors));
	}
	public static double greyscale(BufferedImage OrigImg,int m,int n){
		/**returns the greyscale value of the pixel in the m row and the n column**/
		double[] RGBA = extractRGB(OrigImg,n,m);//n,m are transposed because the input argument of the function is x,y (x for horizontal coordinate and y for vertical)
		double avg = (RGBA[0]+RGBA[1]+RGBA[2])/3;
		return avg;
	}
	public static double oldGreyScaleForForwarding(BufferedImage OrigImg,int m,int n){

		double[] RGBA = extractRGB(OrigImg,m,n);//n,m are not transposed for the use of the costs matrix and forward dynamic maps.
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

