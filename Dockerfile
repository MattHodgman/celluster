FROM python:3
RUN pip3 install pandas
RUN pip3 install numpy

COPY . /app