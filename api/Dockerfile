FROM python:3.7-slim
ADD . /app
WORKDIR /app
RUN pip install -r requirements.txt
EXPOSE 80
CMD quart run --host=0.0.0.0 --port=80