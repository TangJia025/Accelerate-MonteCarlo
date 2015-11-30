public class MonteCarloDemo {    
    public native String nativeMonteCarlo();    

    static {
        System.loadLibrary("montecarlo");
    }        

    public void print () {
    String str = nativeMonteCarlo();
    System.out.println(str);
    }
    
    public static void main(String[] args) {
    (new MonteCarloDemo()).print();
    return;
    }
}
