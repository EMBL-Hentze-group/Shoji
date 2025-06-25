FROM python:3.12 as build
COPY . /shoji
RUN pip install --no-cache-dir poetry \
    && cd shoji \
    && poetry build

FROM python:3.12
COPY --from=build /shoji/dist /tmp/
RUN cd /tmp/ \
    && ls *.whl| xargs -I whl pip install --no-cache-dir whl \
    && rm -rf *.whl \
    && pip cache purge
