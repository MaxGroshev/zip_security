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
import android.provider.Settings
import android.os.Build
import android.database.Cursor
import android.os.Environment
import android.provider.OpenableColumns
import java.io.InputStreamReader
import android.view.View

class MainActivity : AppCompatActivity() {

    private lateinit var tvContent: TextView
    private lateinit var btnOpen: Button
    private lateinit var btnChooseAction: Button
    private lateinit var tvResult: TextView
    private lateinit var tvCurrentlyDisplayingFile: TextView

    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        requestAllFilesAccessPermission()
        setContentView(R.layout.activity_main)

        tvContent = findViewById(R.id.tvChosenFileContent)
        btnOpen = findViewById(R.id.btnChooseFileForArchiving)
        btnChooseAction = findViewById(R.id.btnChooseAction)
        tvResult = findViewById(R.id.tvResult)
        tvCurrentlyDisplayingFile = findViewById(R.id.tvCurrentlyDisplayingFile)

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

    private external fun archiveAndSecure(src: String, dst: String): String
    private external fun unarchiveAndOpen(src: String, dst: String, key: String): String

    enum class ActionType {
        ENCRYPT,
        UNENCRYPT,
        NONE
    }

    private fun showCustomDialog() {
        val dialogView = LayoutInflater.from(this).inflate(R.layout.dialog_action_setup, null)
        val etFileName = dialogView.findViewById<EditText>(R.id.etFileName)
        val etPassword = dialogView.findViewById<EditText>(R.id.etPassword)
        val radioGroup = dialogView.findViewById<RadioGroup>(R.id.radioGroupActions)

        radioGroup.setOnCheckedChangeListener { group, checkedId ->
            when (checkedId) {
                R.id.rbUnencrypt -> {
                    etPassword.visibility = View.VISIBLE
                }
                R.id.rbEncrypt -> {
                    etPassword.visibility = View.GONE
                }
            }
        }
        if (radioGroup.checkedRadioButtonId == R.id.rbEncrypt) {
            etPassword.visibility = View.GONE
        }

        val builder = AlertDialog.Builder(this)
        builder.setView(dialogView)
        builder.setTitle("Параметры операции")

        builder.setPositiveButton("ОК") { dialog, which ->


            val dst = etFileName.text.toString()
            val password = etPassword.text.toString()

            val selectedAction = when (radioGroup.checkedRadioButtonId) {
                R.id.rbEncrypt -> ActionType.ENCRYPT
                R.id.rbUnencrypt -> ActionType.UNENCRYPT
                else -> ActionType.NONE
            }


            if (etFileName.toString().isNotEmpty()) {
                when (selectedAction) {
                    ActionType.ENCRYPT -> {
                        encryptFile(tvCurrentlyDisplayingFile.getText().toString(), dst)
                    }

                    ActionType.UNENCRYPT -> {
                        unencryptFile(tvCurrentlyDisplayingFile.getText().toString(), dst, password)
                    }

                    ActionType.NONE -> {
                        Toast.makeText(this, "Ошибка: действие не выбрано", Toast.LENGTH_SHORT)
                            .show()
                    }
                }
            } else {
                Toast.makeText(this, "Пожалуйста, заполните поля названия файла", Toast.LENGTH_SHORT).show()
            }
        }

        builder.setNegativeButton("Отмена") { dialog, which ->
            dialog.dismiss()
        }

        builder.create().show()
    }

    private fun encryptFile(src: String, dst: String) {
        if (!tvContent.getText().toString().isEmpty()) {
            val download_dir = Environment.getExternalStorageDirectory().getPath() + "/Download/"
            val src_file = download_dir + src;
            val dst_file = download_dir + dst;
            var resultFromCpp = archiveAndSecure(src_file, dst_file)
            tvResult.text = resultFromCpp
        } else {
            tvResult.text = "строка пуста"
        }
    }

    private fun unencryptFile(src: String, dst: String, key: String) {
        if (!tvContent.getText().toString().isEmpty()) {
            val download_dir = Environment.getExternalStorageDirectory().getPath() + "/Download/"
            val src_file = download_dir + src;
            val dst_file = download_dir + dst;
            var resultFromCpp = unarchiveAndOpen(src_file, dst_file, key)
            tvContent.text = resultFromCpp
            tvCurrentlyDisplayingFile.text = dst_file;
        } else {
            tvResult.text = "строка пуста"
        }
    }

    private fun openFile(launcher: androidx.activity.result.ActivityResultLauncher<Intent>) {
        val intent = Intent(Intent.ACTION_OPEN_DOCUMENT).apply {
            addCategory(Intent.CATEGORY_OPENABLE)
            type = "text/plain"
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
            tvCurrentlyDisplayingFile.text = getFileName(uri)
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
                    val nameIndex = it.getColumnIndex(OpenableColumns.DISPLAY_NAME)
                    if (nameIndex >= 0) {
                        result = it.getString(nameIndex)
                    }
                }
            }
        }

        if (result == null) {
            result = uri.path
            val cut = result?.lastIndexOf('/')
            if (cut != null && cut != -1) {
                result = result?.substring(cut + 1)
            }
        }
        return result ?: "unknown_file"
    }

    private fun requestAllFilesAccessPermission() {
        if (Build.VERSION.SDK_INT >= Build.VERSION_CODES.R) {
            if (!Environment.isExternalStorageManager()) {
                try {
                    val intent = Intent(Settings.ACTION_MANAGE_APP_ALL_FILES_ACCESS_PERMISSION)
                    intent.addCategory("android.intent.category.DEFAULT")
                    intent.data = Uri.parse("package:${packageName}")
                    startActivity(intent)
                } catch (e: Exception) {
                    val intent = Intent(Settings.ACTION_MANAGE_ALL_FILES_ACCESS_PERMISSION)
                    startActivity(intent)
                }
            } else {
                Toast.makeText(this, "Access is already granted", Toast.LENGTH_SHORT).show()
            }
        } else {
            throw Error("This line of code should be unreachable")
        }
    }
    companion object {
        init {
            System.loadLibrary("zipper")
        }
    }
}