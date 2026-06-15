@echo off
chcp 65001 >nul
cd /d %~dp0

echo ========================================
echo Plasma cfRNA Brain Injury Tracing System
echo ========================================

python --version >nul 2>nul
if errorlevel 1 (
    echo Python was not detected. Please install Python 3.9 or newer first.
    pause
    exit /b 1
)

python -m pip install -r requirements.txt

echo Opening http://localhost:8501 ...
start "" powershell -NoProfile -WindowStyle Hidden -Command "Start-Sleep -Seconds 4; Start-Process 'http://localhost:8501'"
python -m streamlit run streamlit_app.py --server.headless=true --server.port=8501
pause
