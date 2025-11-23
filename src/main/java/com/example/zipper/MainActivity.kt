package com.example.zipper

import android.app.Activity
import android.content.Intent
import android.net.Uri
import android.os.Bundle
import android.widget.Button
import android.widget.EditText
import android.widget.TextView
import android.widget.Toast
import android.view.LayoutInflater
import android.widget.RadioGroup
import androidx.appcompat.app.AlertDialog
import androidx.appcompat.app.AppCompatActivity
import androidx.activity.result.contract.ActivityResultContracts
import java.io.BufferedReader
import java.io.File
import java.io.IOException
import android.database.Cursor
import android.provider.OpenableColumns
import java.io.InputStreamReader

class MainActivity : AppCompatActivity() {

    private lateinit var tvContent: TextView
    private lateinit var btnOpen: Button
    private lateinit var btnChooseAction: Button
    private lateinit var tvResult: TextView
    private lateinit var tvChosenFileName: TextView

    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        setContentView(R.layout.activity_main)

        tvContent = findViewById(R.id.tvChosenFileContent)
        btnOpen = findViewById(R.id.btnChooseFileForArchiving)
        btnChooseAction = findViewById(R.id.btnChooseAction)
        tvResult = findViewById(R.id.tvResult)
        tvChosenFileName = findViewById(R.id.tvChosenFileName)


        btnChooseAction.setOnClickListener {
            showCustomDialog()
        }

        val openDocLauncher = registerForActivityResult(ActivityResultContracts.StartActivityForResult()) { result ->
            if (result.resultCode == Activity.RESULT_OK) {
                result.data?.data?.let { uri ->
                    readFile(uri)
                }
            }
        }
        btnOpen.setOnClickListener {
            openFile(openDocLauncher)
        }
    }

    private external fun processTextNative(text: String): String
    private external fun archiveAndSecure(text: String, save_to: String): String
    private external fun unarchiveAndOpen(path_to: String): String

    enum class ActionType {
        UPLOAD,
        DOWNLOAD,
        NONE
    }
    private fun showCustomDialog() {
        val dialogView = LayoutInflater.from(this).inflate(R.layout.dialog_action_setup, null)

        val builder = AlertDialog.Builder(this)
        builder.setView(dialogView)
        builder.setTitle("Параметры операции")

        builder.setPositiveButton("ОК") { dialog, which ->
            val etFileName = dialogView.findViewById<EditText>(R.id.etFileName)
            val etPassword = dialogView.findViewById<EditText>(R.id.etPassword)
            val radioGroup = dialogView.findViewById<RadioGroup>(R.id.radioGroupActions)

            val dst = etFileName.text.toString()
            val password = etPassword.text.toString()

            val selectedAction = when (radioGroup.checkedRadioButtonId) {
                R.id.rbUpload -> ActionType.UPLOAD
                R.id.rbDownload -> ActionType.DOWNLOAD
                else -> ActionType.NONE
            }

            when (selectedAction) {
                ActionType.UPLOAD -> {
                    encryptFile(dst, dst, password)
                }
                ActionType.DOWNLOAD -> {
//                    unencryptdFile(src, dst, password)
                }
                ActionType.NONE -> {
                    Toast.makeText(this, "Ошибка: действие не выбрано", Toast.LENGTH_SHORT).show()
                }
            }
        }

        builder.setNegativeButton("Отмена") { dialog, which ->
            dialog.dismiss()
        }

        builder.create().show()
    }

    private fun encryptFile(src: String, dst: String, password: String) {
        if (!tvContent.getText().toString().isEmpty()) {
            val sampleData = """
                    Hellooooo from Android Kotlin and C++!
                    This file was created using JNI.
                    Timestamp: ${System.currentTimeMillis()}
                    Multi-line content works perfectly!
                """.trimIndent()

            val textAsString: String = tvContent.getText().toString()
            val encryped_file = File(filesDir, "encrypte.txt")
            val unencrypted_file = File(filesDir, "unencrypte.txt")

            var resultFromCpp = archiveAndSecure(sampleData, encryped_file.absolutePath)
            resultFromCpp = unarchiveAndOpen(encryped_file.absolutePath)
            tvResult.text = resultFromCpp
        } else {
            tvResult.text = "строка пуста"
        }
    }

    private fun deencryptFile(src: String, dst: String, password: String) {
    }

    // 3. Запуск системного диалога выбора файла
    private fun openFile(launcher: androidx.activity.result.ActivityResultLauncher<Intent>) {
        val intent = Intent(Intent.ACTION_OPEN_DOCUMENT).apply {
            addCategory(Intent.CATEGORY_OPENABLE)
            type = "text/plain" // Фильтр (показывать только текстовые файлы)
        }
        launcher.launch(intent)
    }

    private fun readFile(uri: Uri) {
        try {
            val stringBuilder = StringBuilder()
            contentResolver.openInputStream(uri)?.use { inputStream ->
                BufferedReader(InputStreamReader(inputStream)).use { reader ->
                    var line: String? = reader.readLine()
                    while (line != null) {
                        stringBuilder.append(line).append("\n")
                        line = reader.readLine()
                    }
                }
            }
            tvContent.text = stringBuilder.toString()
            tvChosenFileName.text = getFileName(uri)
        } catch (e: Exception) {
            e.printStackTrace()
            Toast.makeText(this, "Ошибка чтения: ${e.message}", Toast.LENGTH_SHORT).show()
        }
    }

    private fun getFileName(uri: Uri): String {
        var result: String? = null
        if (uri.scheme == "content") {
            val cursor: Cursor? = contentResolver.query(uri, null, null, null, null)
            cursor?.use {
                if (it.moveToFirst()) {
                    // Ищем индекс колонки с именем
                    val nameIndex = it.getColumnIndex(OpenableColumns.DISPLAY_NAME)
                    if (nameIndex >= 0) {
                        result = it.getString(nameIndex)
                    }
                }
            }
        }

        // Если не удалось найти имя через ContentResolver, берем последний сегмент пути
        if (result == null) {
            result = uri.path
            val cut = result?.lastIndexOf('/')
            if (cut != null && cut != -1) {
                result = result?.substring(cut + 1)
            }
        }
        return result ?: "unknown_file"
    }
    companion object {
        // Used to load the 'myapplication' library on application startup.
        init {
            System.loadLibrary("zipper")
        }
    }
}