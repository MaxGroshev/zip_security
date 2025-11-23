package com.example.zipper // Укажите свой пакет

import android.app.Activity
import android.content.Intent
import android.net.Uri
import android.os.Bundle
import android.widget.Button
import android.widget.EditText
import android.widget.TextView
import android.widget.Toast
import androidx.activity.result.contract.ActivityResultContracts
import androidx.appcompat.app.AppCompatActivity
import java.io.BufferedReader
import java.io.File
import java.io.IOException
import java.io.InputStreamReader

class MainActivity : AppCompatActivity() {

    private lateinit var etInput: EditText
    private lateinit var tvContent: TextView
    private lateinit var btnSave: Button
    private lateinit var btnOpen: Button
    private lateinit var btnProcess: Button
    private lateinit var tvResult: TextView

    // Переменная для хранения текста, который мы хотим сохранить
    private var textToSave: String = ""

    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        setContentView(R.layout.activity_main)

        // Инициализация View
        etInput = findViewById(R.id.etInput)
        tvContent = findViewById(R.id.tvContent)
        btnSave = findViewById(R.id.btnSave)
        btnOpen = findViewById(R.id.btnOpen)
        btnProcess =findViewById(R.id.btnProcess)
        tvResult = findViewById(R.id.tvResult)


        // Регистрация "контракта" на создание документа
        val createDocLauncher = registerForActivityResult(ActivityResultContracts.StartActivityForResult()) { result ->
            if (result.resultCode == Activity.RESULT_OK) {
                result.data?.data?.let { uri ->
                    writeFile(uri, textToSave)
                }
            }
        }
        btnSave.setOnClickListener {
            textToSave = etInput.text.toString()
            if (textToSave.isNotEmpty()) {
                createFile(createDocLauncher)
            } else {
                Toast.makeText(this, "Введите текст для сохранения", Toast.LENGTH_SHORT).show()
            }
        }

        // --- ЛОГИКА ОТКРЫТИЯ ---

        // Регистрация "контракта" на открытие документа
        val openDocLauncher = registerForActivityResult(ActivityResultContracts.StartActivityForResult()) { result ->
            if (result.resultCode == Activity.RESULT_OK) {
                result.data?.data?.let { uri ->
                    readFile(uri)
                }
            }
        }

        btnProcess.setOnClickListener {
            if (!tvContent.getText().toString().isEmpty()) {
                val sampleData = """
                    Hellooooo from Android Kotlin and C++!
                    This file was created using JNI.
                    Timestamp: ${System.currentTimeMillis()}
                    Multi-line content works perfectly!
                """.trimIndent()

                val textAsString: String = tvContent.getText().toString()
                val encryped_file = File(filesDir, "encrypted.txt")
                val unencrypted_file = File(filesDir, "unencrypted.txt")


                var resultFromCpp = archiveAndSecure(sampleData, encryped_file.absolutePath)
                resultFromCpp = unarchiveAndOpen(encryped_file.absolutePath)
                tvResult.text = resultFromCpp
            } else {
                tvResult.text = "строка пуста"
            }
        }

        btnOpen.setOnClickListener {
            openFile(openDocLauncher)
        }
    }

    private external fun processTextNative(text: String): String
    private external fun archiveAndSecure(text: String, save_to: String): String
    private external fun unarchiveAndOpen(path_to: String): String


    // Вспомогательная функция для чтения файла из папки assets
    // (Для простоты используем assets, чтобы не возиться с разрешениями Android прямо сейчас)
    private fun readTextFromFile(fileName: String): String {
        return try {
            assets.open(fileName).bufferedReader().use { it.readText() }
        } catch (e: IOException) {
            e.printStackTrace()
            ""
        }
    }
    // 1. Запуск системного диалога создания файла
    private fun createFile(launcher: androidx.activity.result.ActivityResultLauncher<Intent>) {
        val intent = Intent(Intent.ACTION_CREATE_DOCUMENT).apply {
            addCategory(Intent.CATEGORY_OPENABLE)
            type = "text/plain" // MIME тип файла (текст)
            putExtra(Intent.EXTRA_TITLE, "my_note.txt") // Имя файла по умолчанию
        }
        launcher.launch(intent)
    }

    // 2. Запись данных по полученному URI
    private fun writeFile(uri: Uri, content: String) {
        try {
            contentResolver.openOutputStream(uri)?.use { outputStream ->
                outputStream.write(content.toByteArray())
            }
            Toast.makeText(this, "Файл успешно сохранен!", Toast.LENGTH_SHORT).show()
        } catch (e: Exception) {
            e.printStackTrace()
            Toast.makeText(this, "Ошибка сохранения: ${e.message}", Toast.LENGTH_SHORT).show()
        }
    }

    // 3. Запуск системного диалога выбора файла
    private fun openFile(launcher: androidx.activity.result.ActivityResultLauncher<Intent>) {
        val intent = Intent(Intent.ACTION_OPEN_DOCUMENT).apply {
            addCategory(Intent.CATEGORY_OPENABLE)
            type = "text/plain" // Фильтр (показывать только текстовые файлы)
        }
        launcher.launch(intent)
    }

    // 4. Чтение данных по URI
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
        } catch (e: Exception) {
            e.printStackTrace()
            Toast.makeText(this, "Ошибка чтения: ${e.message}", Toast.LENGTH_SHORT).show()
        }
    }
    companion object {
        // Used to load the 'myapplication' library on application startup.
        init {
            System.loadLibrary("zipper")
        }
    }
}