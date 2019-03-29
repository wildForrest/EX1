
import java.io.File;
import java.io.IOException;
import java.awt.image.BufferedImage;

import javax.imageio.ImageIO;

public class seamCarving {
	
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
		 File f=new File("images\\lake.jpg");
			BufferedImage img=null;
			
			try {
				 img = ImageIO.read(f);
			} catch (IOException e) {
				
				e.printStackTrace();
			}
			
			int width = img.getWidth();
			int height = img.getHeight();
			double energy_map[][] = createEnegryMap(img, width,height);
			/*System.out.println("img height = "+img.getHeight()+
					" img width ="+img.getWidth());
			System.out.println("energy height = "+energy_map.length+
					" energy width ="+energy_map[0].length);*/
			
			ExportMatrixAsImage(energy_map,img,"images\\energymap.jpg");
				
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

	public static double[][] addEntropy(double[][] Emap,BufferedImage img){
		/* part 1 - 2
		 * roee*/
		
		return null;
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
