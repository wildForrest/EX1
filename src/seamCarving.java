
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.awt.image.BufferedImage;

import javax.imageio.ImageIO;

public class seamCarving {
	public static int Pcounter = 0;//number of calculation spared by memoization of p-values
	public static int Gcounter = 0;//number of calculation spared by memoization of greyscale
	
	public static void main(String args[]) {
		tmpTest();
	}

	public static void exportEmptyPicture(){
		try{
		    BufferedImage image = null;
			image = new BufferedImage(940, 940, BufferedImage.TYPE_INT_ARGB);
		      File f = new File("images\\Empty.jpg");  //output file path
		      ImageIO.write(image, "jpg", f);
		      System.out.println("Writing complete.");
		    }catch(IOException e){
		      System.out.println("Error: "+e);
		    }
	}
	public static void tmpTest(){
		 File f=new File("images\\strawberry.jpg");
			BufferedImage img=null;
			
			try {
				 img = ImageIO.read(f);
			} catch (IOException e) {
				
				e.printStackTrace();
			}
			
			int width = img.getWidth();
			int height = img.getHeight();
			double energy_map[][] = createEnegryMap(img, width,height);
			System.out.println("finish creating energy map");
			
			ExportMatrixAsImage(energy_map,img,"images\\energyMap.jpg");
			double energy_map_after_entropy[][] =addEntropy(energy_map,img);
			ExportMatrixAsImage(energy_map_after_entropy,img,"images\\energyMapWithEntropy.jpg");
			System.out.println("finish adding entropy");
			System.out.println(Pcounter);
			System.out.println(Gcounter);	
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
		/**author: Roee
		 * the indices i,j are given in the opposite order in the arguments list
		 * because setRGB(i, j, p) approaches the TRANSPOSE matrix of the picture**/
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

	public static double[][] addEntropy(double[][] Emap,BufferedImage OrigImg){
		/**author: Roee
		 * part 1 - 2**/
		HashMap<Integer,HashMap<Integer,Double>> pValuesDictionary = new HashMap<>();
		HashMap<Integer,HashMap<Integer,Double>> greyscaleDictionary = new HashMap<>();
		int rows = Emap.length;
		int columns = Emap[0].length;
		double[][] res = new double[rows][columns];
		for(int i=0;i<rows;i++){
			for(int j=0;j<columns;j++){
				res[i][j]=Emap[i][j]+pixelEntropy(OrigImg,i,j,
						pValuesDictionary,greyscaleDictionary);
			}
		}
		
		
		return res;
	}
	public static double pixelEntropy(BufferedImage OrigImg,int i,int j,
			HashMap<Integer,HashMap<Integer,Double>> pValuesDictionary,HashMap<Integer,HashMap<Integer,Double>> greyscaleDictionary){
		int rows = OrigImg.getHeight();
		int columns = OrigImg.getWidth();

		int left = Math.max(0, i-4);
		int right = Math.min(columns-1, i+4);
		int up = Math.max(0, j-4);
		int down = Math.min(rows-1, j+4);
		double sum=0;
		double p;
		for(int m=left;m<right;m++){
			for(int n=up;n<down;n++){
				if(dictionaryContainsKey(pValuesDictionary,m,n)){
					p= dictionaryGet(pValuesDictionary,m,n);
					Pcounter++;
				}
				else{
					p=p(OrigImg,m,n,greyscaleDictionary);
					dictionaryPut(pValuesDictionary,m,n,p);
					}
				sum+=p*Math.log(p);
			}
		}
		return -sum;
	}
	public static double p(BufferedImage OrigImg,int m,int n,
			HashMap<Integer,
			HashMap<Integer,Double>> greyscaleDictionary){
		int rows = OrigImg.getHeight();
		int columns = OrigImg.getWidth();
		int left = Math.max(0, m-4);
		int right = Math.min(columns-1, m+4);
		int up = Math.max(0, n-4);
		int down = Math.min(rows-1, n+4);
		double greyscale;
		double sum=0;
		for(int k=left;m<right;m++){
			for(int l=up;n<down;n++){
				if(dictionaryContainsKey(greyscaleDictionary,k,l)){
					greyscale= dictionaryGet(greyscaleDictionary,k,l);
					Gcounter++;
				}
				else{
					greyscale=greyscale(OrigImg,k,l);
					dictionaryPut(greyscaleDictionary,k,l,greyscale);
			}
				sum+=greyscale;
		}}
		
		return greyscale(OrigImg,m,n)/sum;
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
	public static double[][] dynamicMap( double[][] Emap){
		/* part 1 - 3
		 * rachel*/
		return null;
	}
	public static int[] chooseSeam(double[][] Dmap, int direction){
		/* part 1 - 4
		 * direction = 0 for horizontal, 1 for verticle
		 * 
		 * */
		return null;
	}
	public static BufferedImage deleteSeam(BufferedImage img, int direction, int[] seamVector){
		/* part 1 - 4
		 * direction = 0 for horizontal, 1 for verticle
		 * 
		 * */
		return null;
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
