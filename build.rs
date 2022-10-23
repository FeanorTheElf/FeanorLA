
fn main() {
    println!("cargo:rustc-link-search=static=mpir");
    println!("cargo:rustc-link-search=native=E:/Users/Simon/Documents/Projekte/mpir/msvc/vs19/lib_mpir_gc/x64/Release");
}