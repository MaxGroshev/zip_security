package nativecpp

class LzwArchiver {
    public external fun archiveAndSecure(src: String, dst: String): String
    public external fun unarchiveAndOpen(src: String, dst: String, key: String): String

    companion object {
        init {
            System.loadLibrary("zipper")
        }
    }
}