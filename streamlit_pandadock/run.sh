#!/bin/bash

# PandaDock Streamlit Application Runner
# This script starts the Streamlit web interface for PandaDock

echo "🐼 Starting PandaDock Streamlit Web Interface..."
echo "=================================================="

# Check if streamlit is installed
if ! command -v streamlit &> /dev/null; then
    echo "❌ Streamlit is not installed. Installing now..."
    pip install streamlit
fi

# Check if required packages are installed
echo "📦 Checking dependencies..."
pip install -r requirements.txt

# Start the Streamlit application
echo "🚀 Launching PandaDock Web Interface..."
echo "📱 The application will open in your default web browser"
echo "🌐 URL: http://localhost:8501"
echo "⏹️  Press Ctrl+C to stop the application"
echo "=================================================="

streamlit run app.py