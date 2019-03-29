import java.io.File;
import java.io.IOException;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;

public class graphicTests {
	public static void main(String args[]){
		readAndWriteAnImage();
	   
	  

}
	 public static void readAndWriteAnImage(){
		 int width = 340;    //width of the image
		    int height = 940;   //height of the image
		    BufferedImage image = null;
		    File f = null;

		    //read image
		    try{
		      f = new File("strawberry.jpg"); //image file path
		      image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		      image = ImageIO.read(f);
		      System.out.println("Reading complete.");
		    }catch(IOException e){
		      System.out.println("Error: "+e);
		    }
		    
		  //write image
		    try{
		      f = new File("strawberry2.jpg");  //output file path
		      ImageIO.write(image, "jpg", f);
		      System.out.println("Writing complete.");
		    }catch(IOException e){
		      System.out.println("Error: "+e);
		    }
		    
	    }
	
}
