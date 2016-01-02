.onLoad = function(libname = find.package("taoR"), pkgname = "taoR") {
    tao_init();
}

.onUnload = function(libpath) {
    tao_finalize();
}