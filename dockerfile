FROM python:3.12.9-slim

ENV PACKAGES "pandas numpy matplotlib seaborn openpyxl"

RUN pip install --upgrade pip

RUN pip install $PACKAGES