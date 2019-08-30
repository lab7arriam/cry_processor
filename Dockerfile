FROM python:3

RUN apt update && apt install -y make cmake gcc build-essential && \
	pip install biopython && \
	rm -rf /var/lib/apt/lists/* && \
	mkdir /cry_processor

ENV PATH=$PATH:/cry_processor

WORKDIR /data

ADD . / /cry_processor/

RUN chmod 775 /cry_processor/cry_processor.py

CMD cry_processor.py
